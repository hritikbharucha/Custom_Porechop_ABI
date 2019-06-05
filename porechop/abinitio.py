"""
Ab-Initio wrapper created by Quentin Bonenfant(quentin.bonenfant@gmail.com
                                       quentin.bonenfant@univ-lille.fr)
https://github.com/qbonenfant

This script is a "wrapper" to use my algorithm to find adapter sequence from the reads instead of
using the adapter.py static database. This is still in its early stage, but do seems to work fairly well.

A futur, far cleaner version of this script will be edited in order to use my c++ code "the proper way".

In order to change the parameters for adaptFinder, change the abinitio.conf file.
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


def greedy_assembl(g):
    """Greedy assembly method to compute the adapter
    TODO: modify it so it use directly the graph...
    @param the overlap graph of kmers
    @return the longest debruijn sequence starting by the first kmer
    """
    kmer_dict = g.nodes(data=True)
    kmer_list = list( dict(kmer_dict).keys() )
    kmer_list.sort( key= lambda x: kmer_dict[x]['weight'])
    kmer_list.reverse()

    km = kmer_list[0]
    ov = km
    found = True
    used = [km]
    #annotating graph
    g.node[km]["path"] = (g.node[km]["path"] + ",greedy").lstrip(",")

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
                    g.node[km2]["path"] = (g.node[km2]["path"] + ",greedy").lstrip(",")
                    break

                elif(reverse != "" and direct == ""):   
                    ov = reverse
                    found = True
                    used.append(km2)
                    g.node[km2]["path"] = (g.node[km2]["path"] + ",greedy").lstrip(",")

                    break
    return(ov)


def dag_heaviest_path(G):
    """Returns the heaviest path in a DAG

    Parameters
    ----------
    G : NetworkX DiGraph
        Graph

    Returns
    -------
    path : list
        Heaviest path

    Comment
    -------
    This is a modified version of the dag_longest_path
    using node weight as distance.
    """
    dist = {}  # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0] + G.node[v]["weight"], v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node, (length, _) = max(dist.items(), key=lambda x: x[1])
    path = []
    while length > 0:
        path.append(node)
        length, node = dist[node]
    return list(reversed(path))


def heavy_path(g):
    """ Searching the truly heaviest path between all source and target nodes.
        Even if longer path tend to be heavier, very heavy short path can 
        also be selected.
    """
    
    hv_path = dag_heaviest_path(g)

    for n in hv_path:
        g.node[n]["path"] = (g.node[n]["path"] + ",heavy").lstrip(",")

    return( hv_path[0][:-1] + "".join( el[-1] for el in hv_path ) )


def longest_path(g):
    lg_path = nx.dag_longest_path(g)

    for n in lg_path:
        g.node[n]["path"] = (g.node[n]["path"] + ",long").lstrip(",")

    return(lg_path[0][:-1] + "".join( el[-1] for el in lg_path ))


def build_graph(count_file):
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
            g.add_node(km, weight = int(nb))#, path = "" )
        
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

    # Returning only the biggest connected component
    return( g.subgraph(max(nx.weakly_connected_components(g), key= lambda x: get_weight(g,x))) )


def execFindAdapt(args):
    
    """ Launch adaptFinder to perform the approximate kmers count and try 
    to rebuild the adapter using different methods.
        The count will be stored as a temporary file in a ./tmp folder.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    # getting args
    fasta_file = args.input
    just_print = args.guess_adapter_only
    print_dest = args.print_dest
    v = args.verbosity


    greedy_adapter = {}
    long_adapter  = {}
    heavy_adapter  = {}

    # Temporary files prefix
    filename_pref = "./tmp/temp_file_"
    tmpDir =  os.path.dirname(filename_pref) 
    
    # Creating tmp folder if not existing
    if( not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)

    # Searchng path to adaptFinder
    adapt_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "adaptFinder")

    out_file_name =  filename_pref + "approx_kmer_count"
    nb_thread = str(min( 4, cpu_count() ))

    command =  adapt_path +" "+ fasta_file + " -k 16 -lc 1.2 -sn 20000 -sl 80 --limit 1000 -nt " + nb_thread + " -v " + str(v) + " -o " + out_file_name
    
    try:    
        subprocess.check_call( command.split() )
    except SystemError as e:
        print(e)
        sys.exit("ERROR: approximate k-mer count failed")

    # Building adapters from counts using  different method:
    # - Building a Deruijn graph and searching the longest path and heaviest path
    # - Greedy assembly based on kmer rank (most frequent first)
    # TODO: I REALLY need to improve the reconstruction using proper assembly method...
    # The proposed methods slightly overstimate the real adapter length and may prefer 
    # adapter containing insertion errors.
    # Porechop can still properly detect adapters since it uses a quite low identity threshold of 75%. 


    for which_end in ["start","end"]:
        if(v>=1):
            print("Assembling " + which_end + " adapters", file=print_dest)
        # Building graph
        g = build_graph(out_file_name + "." + which_end )

        # preping for anotation
        nx.set_node_attributes(g, "", "path" )

        # adapters
        if(v>=1):
            print("\tBuilding greedy adapter", file = print_dest)
        greedy_adapter[which_end] = greedy_assembl(g)
        
        try:
            if(v>=1):
                print("\tBuilding longest path adapter", file = print_dest)
            long_adapter[which_end] = longest_path(g)

        # A lot of exception could append here...
        except:
            print("\tCould not compute " + which_end + " adaper using longest path  method", file = print_dest)
            print("\tThe resulting graph probably contains a loop.", file = print_dest)
            long_adapter[which_end]  = ""
            heavy_adapter[which_end] = ""

        if(v>=1):
            print("\tBuilding heavy path adapter", file = print_dest)            
        heavy_adapter[which_end]  = heavy_path(g)

        #Exporting, if required
        if( args.export_graph is not None):
            if(v>=1):
                print("\tExporting assembly graph", file = print_dest)
            
            base = "_adapter_graph"
            if( ".graphml" in args.export_graph):
                base = "_" + os.path.basename(args.export_graph)
            else:
                path = os.path.join( os.path.dirname(args.export_graph), which_end + base + ".graphml") 
            
            nx.write_graphml(g,path)




    # If we just need to print the adapter
    if(just_print):
        if(v>0):
            print("\n\nINFERRED ADAPTERS:\n", file= print_dest )

            print("Greedy assembly method", file= print_dest )
            print("Start:\t" + greedy_adapter["start"], file= print_dest )
            print("End:\t"   + greedy_adapter["end"]+"\n", file= print_dest )

            print("Longest Path assembly method", file= print_dest )
            print("Start:\t" + long_adapter["start"], file= print_dest )
            print("End:\t"   + long_adapter["end"]+"\n", file= print_dest )

            print("Heaviest Path assembly method", file= print_dest )
            print("Start:\t" + heavy_adapter["start"], file= print_dest )
            print("End:\t"   + heavy_adapter["end"]+"\n", file= print_dest )
        # if non verbose mode, display should be minimal
        else:
            print("greedy")
            print(greedy_adapter["start"])
            print(greedy_adapter["end"])
            print("longest_path")
            print(long_adapter["start"])
            print(long_adapter["end"])
            print("heaviest_path")
            print(heavy_adapter["start"])
            print(heavy_adapter["end"])


    # If we need to use them, getting back to porechop objects
    else:
        if(v>=1):
            print("Building adapter object", file= print_dest )
        adp = [Adapter("abinitio_long_adapter",
            start_sequence=('abinitio_graph_Top' , long_adapter['start']),
            end_sequence=('abinitio_graph_Bottom', long_adapter['end'])),

            Adapter("abinitio_heavy_adapter",
            start_sequence=('abinitio_heavy_Top' , heavy_adapter['start']),
            end_sequence=('abinitio_heavy_Bottom', heavy_adapter['end'])),

            Adapter("abinitio_greedy_adapter",
            start_sequence=('abinitio_greedy_Top' , greedy_adapter['start']),
            end_sequence=('abinitio_greedy_Bottom', greedy_adapter['end']))
        ]
        if(v>=1):
            print("The inference of adapters sequence is done.", file= print_dest )
        return(adp)
    return([])


if __name__ == '__main__':
    adapt = execFindAdapt(sys.argv[1])
    for ad in adapt:
        print(ad.__dict__)