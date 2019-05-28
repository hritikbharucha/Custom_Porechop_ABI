"""
Ab-Initio wrapper created by Quentin Bonenfant(quentin.bonenfant@gmail.com
                                       quentin.bonenfant@univ-lille.fr)
https://github.com/qbonenfant

This script is a "wrapper" to use my algorithm to find adapter sequence from the reads instead of
using the adapter.py static database. This is still in its early stage, but do seems to work fairly well.

A futur, far cleaner version of this script will be edited in order to use my c++ code "the proper way".

I could have allowed users to specify arguments like k-mer size, LC threshold, sample size and length,
but i will need to alter the porechop parser a lot more.
"""

import os
import subprocess
from multiprocessing import cpu_count
from  .adapters import Adapter
import networkx as nx
import sys

def haveOverlap(seq1, seq2):
    """Check if the sequence 1 is a prefix of sequence 2
    @param first sequence
    @param second sequence
    @return is there an overlap ?
    """
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")

def get_weight(g,path):
    """Compute the weight of a path in the graph
    @param the graph
    @param the path (list of nodes)
    @return the total weight of this path

    """
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)

def greedy_assembl(kmer_list):
    """Greedy assembly method to compute the adapter
    @param the list of kmer, ordered by count (descending)
    @return the longest debruijn sequence starting by the first kmer
    """

    km = kmer_list[0]
    ov = km
    found = True
    used = [km]
    while(found):
        found = False
        for km2 in kmer_list:
            if(km2 not in used):
                direct = haveOverlap(ov, km2)
                reverse = haveOverlap(km2, ov)
                if(direct != "" and reverse == ""):
                    ov = direct
                    found = True
                    used.append(km2)
                    break

                elif(reverse != "" and direct == ""):   
                    ov = reverse
                    found = True
                    used.append(km2)
                    break
    return(ov)



def heavy_path(g):
    """ Searching the truly heaviest path between all source and target nodes.
        Even if longer path tend to be heavier, very heavy short path can 
        also be selected.
    """

    sources = [n for n in g.nodes if not list(g.predecessors(n))  ]
    targets = [n for n in g.nodes if not list(g.successors(n))    ]
    # print(sources)
    # print(targets)
    hv_path = []
    w = 0
    for source in sources:
        for target in targets:
            try:
                heaviest_path = max((path for path in nx.all_simple_paths(g, source, target)),
                            key=lambda path: get_weight(g,path))
                current_w = get_weight(g,heaviest_path)
                if(current_w > w):
                    hv_path = heaviest_path
                    w = current_w
            # I know it's bad, i will fix this.
            #TODO: fix this
            except ValueError:
                pass
    return( hv_path[0][:-1] + "".join( el[-1] for el in hv_path ) )




def buildAdapter(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
    """

    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # building graph
    with open(count_file, 'r') as f:
        for line in f:
            km, nb = line.rstrip("\n").split("\t")
            kmer_count[km] = int(nb)
            kmer_list.append(km)
            g.add_node(km, weight = int(nb), path = "" )
        
    # searching overlaps
    for km in kmer_list:
        for km2 in kmer_list:
            # since i do all the possible combinations
            # There is no need to test both orientation
            direct = haveOverlap(km, km2)
            if(direct != ""):            
                g.add_edge(km, km2)
            
    # removing singletons
    g.remove_nodes_from( [node for node in g.nodes if  len(list(nx.all_neighbors(g,node))) == 0] )
    
    # Longest path        
    lg_path = nx.dag_longest_path(g)

    # adapters
    graph_adapter  = lg_path[0][:-1] + "".join( el[-1] for el in lg_path )
    greedy_adapter = greedy_assembl(kmer_list)
    heavy_adapter  = heavy_path(g)

    return( graph_adapter, greedy_adapter, heavy_adapter )


def execFindAdapt(fasta_file):
    
    """ Launch adaptFinder to perform the approximate kmers count and try 
    to rebuild the adapter using different methods.
        The count will be stored as a temporary file in a ./tmp folder.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    greedy_adapter = {}
    graph_adapter  = {}
    heavy_adapter  = {}

    # Temporary files prefix
    filename_pref = "./tmp/temp_file_"
    tmpDir =  os.path.dirname(filename_pref) 
    
    # Creating tmp folder if not existing
    if( not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)


    # Working on both ends
    for which_end in ['starts','ends']:

        #Counting k-mers at ~ 2 errs
        print("Counting kmers at the " + which_end + " of the sequences")

        adapt_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "adaptFinder")

        bot = "" if which_end == "starts" else " -bot " 
        command =  adapt_path +" "+ fasta_file + " -k 16 -lc 1.1 -sn 30000 -sl 80 --limit 1000 -nt " + str(min( 4, cpu_count() )) + bot + " -o " + filename_pref + which_end + "_count.txt"
        subprocess.check_call( command.split() )

        # Building adapters from counts using  different method:
        # - Building a Deruijn graph and searching the longest path and heaviest path
        # - Greedy assembly based on kmer rank
        # TODO: I REALLY need to improve the reconstruction using proper assembly method...
        # The proposed methods slightly overstimate the real adapter length and may prefer 
        # adapter containing insertion errors.
        # Porechop can still properly detect adapters since it uses a quite low identity threshold of 75%. 
        
        graph, greed, heavy = buildAdapter(filename_pref + which_end + "_count.txt")
        
        graph_adapter[which_end]  = graph
        greedy_adapter[which_end] = greed
        heavy_adapter[which_end]  = heavy

    # Getting back to porechop objects
    print("Building adapter object")
    adp = [Adapter("abinitio_graph_adapter",
        start_sequence=('abinitio_graph_Top' , graph_adapter['starts']),
        end_sequence=('abinitio_graph_Bottom', graph_adapter['ends'])),

        Adapter("abinitio_heavy_adapter",
        start_sequence=('abinitio_heavy_Top' , heavy_adapter['starts']),
        end_sequence=('abinitio_heavy_Bottom', heavy_adapter['ends'])),

        Adapter("abinitio_greedy_adapter",
        start_sequence=('abinitio_greedy_Top' , greedy_adapter['starts']),
        end_sequence=('abinitio_greedy_Bottom', greedy_adapter['ends']))
    ]

    print("The inference of adapters sequence is done.")
    
    return(adp)


if __name__ == '__main__':
    adapt = execFindAdapt(sys.argv[1])
    for ad in adapt:
        print(ad.__dict__)