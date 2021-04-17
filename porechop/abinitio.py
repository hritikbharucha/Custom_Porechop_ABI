"""
Ab-Initio wrapper created by Quentin Bonenfant
quentin.bonenfant@univ-lille.fr
quentin.bonenfant@gmail.com

https://github.com/qbonenfant
https://github.com/bonsai-team

This script is a "wrapper" to use my algorithm to find adapter sequence
from the reads instead of using the adapter.py static database. This is
still in its early stage, but do seems to work fairly well.
A futur, far cleaner version of this script will be edited in order to use
my c++ code "the proper way".
In order to change the parameters for adaptFinder,
change the porchop/abinitio.conf file.

This works has been funded by the french National Ressearch Agency (ANR).
"""

import os
import subprocess
from multiprocessing import cpu_count
from .adapters import Adapter
from collections import defaultdict as dd
from statistics import median, mean
import networkx as nx
import sys


##############################################################################
#                                 CONSTANTS                                  #
##############################################################################

CUT_RATIO = 0.075  # /!\ The way this value is has been changed.
METHODS = ["greedy", "heavy"]


##############################################################################
#                              UTILITY FUNCTIONS                             #
##############################################################################


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


def get_weight(g, path):
    """Compute the weight of a path in the graph
    @param the graph
    @param the path (list of nodes)
    @return the total weight of this path

    """
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)


def concat_path(path):
    """Concat the kmers of a path into a single sequence
    @param The path as a k-mer list
    @return the full sequence.
    """
    return(path[0][:-1] + "".join(el[-1] for el in path))


def check_drop(path, g):
    """DEPRECATED - see start_cut()
    Check if there is a frequency drop in the path
    and cut the adapter if necessary.
    This function cut the end of forward adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @return The adjusted adapter path
    """
    last = 0
    cut = 0
    for i in range(len(path) - 1, 1, -1):
        last = path[i]
        prev = path[i - 1]
        if(g.nodes[prev]["weight"] / float(g.nodes[last]["weight"]) > CUT_RATIO):
            cut -= 1
        else:
            break
    if(cut == 0):
        return(path)
    return(path[:cut])


def check_drop_back(path, g):
    """DEPRECATED - see end_cut()
    Check if there is a frequency drop in the path
    and cut the adapter if necessary.
    This function cut the start of reverse adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @return The adjusted adapter path
    """
    first = 0
    cut = 0
    for i in range(len(path) - 1):
        first = path[i]
        nxt = path[i + 1]
        if(g.nodes[nxt]["weight"] / float(g.nodes[first]["weight"]) > CUT_RATIO):
            cut += 1
        else:
            break
    return(path[cut:])


def start_cut(start_adp_data, g, w=7):
    """ Try to find the best place to cut the raw adapter path
     to get an appropriate adapter.
    This function cut the end of forward adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @param (opt) w, window size for sliding median smoothing (default: 7)
    @return The adjusted adapter path
    """

    # Data for Start Adapters
    start_dist = [g.nodes[el]["weight"] for el in start_adp_data]

    # computing start epsilon based on counts
    start_epsilon = CUT_RATIO * max(start_dist)

    # Sliding median method
    start_sl_md = sliding_median(start_dist, w)
    # # Sliding mean method
    # start_sl_md = sliding_mean(start_dist, w)

    # Simple derivative
    start_sl_md = simple_deriv(start_sl_md)

    # Finding cutting point
    # start median zone method
    smz = find_median_zone(start_sl_md, start_epsilon)

    # Finding the position to cut
    start_cut = len(start_adp_data) - 1
    if(smz):
        start_cut = min(max(smz), start_cut)

    # Returning the cut path
    return(start_adp_data[0:start_cut])


def end_cut(end_adp_data, g, w=7):
    """ Try to find the best place to cut the raw adapter path
     to get an appropriate adapter.
    This function cut the end of end adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @param (opt) w, window size for sliding median smoothing (default: 7)
    @return The adjusted adapter path
    """

    # Fetching weight distribution
    end_dist = [g.nodes[el]["weight"] for el in end_adp_data]

    # Computing end epsilon based on counts
    end_epsilon = CUT_RATIO * max(end_dist)

    # Sliding median method
    end_sl_md = sliding_median(end_dist, w)
    # # Sliding mean method
    # end_sl_md = sliding_mean(end_dist, w)

    # Simple derivative
    end_sl_md = simple_deriv(list(reversed(end_sl_md)))

    # Finding cutting point:
    # end median zone method
    emz = find_median_zone(end_sl_md, end_epsilon)

    # re-reversing the end array
    end_sl_md.reverse()
    emz = [len(end_sl_md) - (s + 1) for s in emz]

    # Finding the position to cut
    end_cut = 0
    if(emz):
        end_cut = max(0, min(emz) + 1)

    # Returning the cut path
    return(end_adp_data[end_cut:])


def sliding_median(counts, w=7):
    lc = len(counts)
    offset = w // 2  # offset to work on "middle" position
    results = []
    for i in range(lc):
        start = max(0, i - offset)
        end = min(i + offset, lc)
        # getting on middle position
        med_rng = counts[start: end]
        results.append(median(med_rng))
    return(results)


def sliding_mean(counts, w=7):
    lc = len(counts)
    offset = w // 2  # offset to work on "middle" position
    results = []
    for i in range(lc):
        start = max(0, i - offset)
        end = min(i + offset, lc)
        # getting on middle position
        med_rng = counts[start: end]
        results.append(mean(med_rng))
    return(results)


def simple_deriv(counts):
    lc = len(counts)
    results = [0]
    for i in range(1, lc, 1):
        val = counts[i] - counts[i - 1]
        results.append(val)
    return(results)


def find_median_zone(counts, epsilon):
    c_median = median(counts)
    s_rng = 0  # start of low zone
    is_zone = False
    s_zone = []
    for i in range(len(counts) - 1, -1, -1):
        c = counts[i]
        if(c < c_median - epsilon):
            is_zone = True
            s_rng = i
                
        else:
            if(is_zone):
                is_zone = False
                s_zone.append(s_rng)
    return(s_zone)


def print_result(adapters, v=1, print_dest=sys.stdout):
    """Display the result, adjusting format depending on selected verbosity
    @param an adapter dictionnary, using appropriate methods as keys
    @param v, the verbosity level, default to 1
    @param print destination, can be custom, but stdout is default
    """

    global METHODS

    out = print_dest

    if(v > 0):
        print("\n\nINFERRED ADAPTERS:\n",
              file=out)
    else:
        out = sys.stdout

    for meth in METHODS:
        msg = meth
        srt = adapters[meth]["start"]
        end = adapters[meth]["end"]
        if(srt and end):
            if(v >= 1):
                meth += " assembly method"
                meth = meth.capitalize()
                srt = "Start:\t" + srt
                end = "End:\t" + end
            print(msg, file=out)
            print(srt, file=out)
            print(end, file=out)


##############################################################################
#                                 ASSEMBLY                                   #
##############################################################################

def greedy_assembl(g):
    """Greedy assembly method to compute the adapter using the graph as input.
    @param the De Bruijn graph of kmers
    @return the longest debruijn sequence starting by the first kmer
    """
    start = max(g.nodes, key=lambda x: g.nodes[x]["weight"])
    path = [start]

    right_node = start
    left_node = start
    while(left_node or right_node):

        # forward extension
        if(right_node):
            r_list = list(g.successors(right_node))
            if(not r_list):
                right_node = None
            else:
                right_node = max(r_list, key=lambda x: g.nodes[x]["weight"])
                path.append(right_node)

        # reverse extension
        if(left_node):
            l_list = list(g.predecessors(left_node))

            if(not l_list):
                left_node = None
            else:

                left_node = max(l_list, key=lambda x: g.nodes[x]["weight"])
                path = [left_node] + path
    return(path)


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
        pairs = [(dist[v][0] + G.nodes[v]["weight"], v) for v in G.pred[node]]
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

    # annotating graph
    for n in hv_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",heavy").lstrip(",")

    return(hv_path)


def longest_path(g):
    """ DEPRECATED - This function is not used anymore
        Searching the longest path between all source and target nodes.
    """
    lg_path = nx.dag_longest_path(g)

    # annotating graph
    for n in lg_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",long").lstrip(",")

    return(lg_path)


def greedy_path(g):
    """Greedy assembly of the adapter from the graph,
    starting from the most frequent
    """

    gd_path = greedy_assembl(g)
    # annotating graph
    for n in gd_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",greedy").lstrip(",")
    return(gd_path)

##############################################################################
#                              GRAPH BUILDING                                #
##############################################################################


def build_graph(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
    """

    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # Avoiding IO error
    try:
        # building graph
        with open(count_file, 'r') as f:
            for line in f:
                km, nb = line.rstrip("\n").split("\t")
                kmer_count[km] = int(nb)
                kmer_list.append(km)
                g.add_node(km, weight=int(nb))  # , path = "" )

    except FileNotFoundError:
        print("\n/!\\ Unable to open k-mer count file:", file=sys.stderr)
        print(count_file, file=sys.stderr)
        print("Either the file was moved, deleted, or filename is invalid.",
              file=sys.stderr)
        print("It is also possible that end adapter ressearch was skipped.",
              file=sys.stderr)
        print("Be sure skip_end / se option is deactivated in adaptFinder.\n",
              file=sys.stderr)

    else:
        # searching overlaps
        for km in kmer_list:
            for km2 in kmer_list:
                # since i do all the possible combinations
                # There is no need to test both orientation
                direct = haveOverlap(km, km2)
                if(direct != ""):
                    g.add_edge(km, km2)

        # removing singletons
        g.remove_nodes_from([node for node in g.nodes if len(
            list(nx.all_neighbors(g, node))) == 0])

        # Returning only the biggest connected component
        g = g.subgraph(max(nx.weakly_connected_components(g),
                           key=lambda x: get_weight(g, x)))
    finally:
        return(g)

##############################################################################
#                              ADPATER BUILDING                              #
##############################################################################


def build_adapter(out_file_name, args):
    """Building adapters from counts using different method:
        - Building a Debruijn graph and searching heaviest path
        - Greedy assembly based on kmer rank (most frequent first)
    Porechop can still properly detect adapters if there is a slight
    difference  since it uses a quite low identity threshold of 75%.
    @param A De Bruijn graph with kmer count as weight
    @param The argparse argument dictionnary from porechop.py
    @return A list of Adapter objects containing infered adapters.
    """

    just_print = args.guess_adapter_only
    print_dest = args.print_dest
    v = args.verbosity

    adapters = dd(dict)  # temporary adapter string dict
    adp = []  # final adpter object list

    # Failsafe if the graphs for the adapter can not be built
    unable_to_build = False

    for which_end in ["start", "end"]:
        if(v >= 1):
            print("Assembling " + which_end + " adapters", file=print_dest)

        # Building graph
        g = build_graph(out_file_name + "." + which_end)

        # Checking graph is properly build
        if(not nx.is_empty(g)):

            # preping for anotation
            nx.set_node_attributes(g, "", "path")

            # Greedy Adapter
            if(v >= 1):
                print("\tBuilding greedy adapter", file=print_dest)

            greedy_p = greedy_path(g)
            cut_greedy_p = []
            if(which_end == "start"):
                cut_greedy_p = start_cut(greedy_p, g)
            elif(which_end == "end"):
                cut_greedy_p = end_cut(greedy_p, g)
            # If the greedy path is not empty, concat it in a string
            if(cut_greedy_p):
                adapters["greedy"][which_end] = concat_path(cut_greedy_p)
            else:
                adapters["greedy"][which_end] = ""

            # Heavy adapter
            if(v >= 1):
                print("\tBuilding heavy path adapter", file=print_dest)
            try:
                heavy_p = heavy_path(g)
                cut_heavy_p = []
                if(which_end == "start"):
                    cut_heavy_p = start_cut(heavy_p, g)
                elif(which_end == "end"):
                    cut_heavy_p = end_cut(heavy_p, g)

                # If the heavy path is not empty, concat it in a string
                if(cut_heavy_p):
                    adapters["heavy"][which_end] = concat_path(cut_heavy_p)
                else:
                    adapters["heavy"][which_end] = ""

            # I know i should specify the exception, but nx exception seems
            # to not be caught if specified here...
            except (nx.NetworkXNotImplemented, nx.HasACycle) as nxE:
                print("\t/!\\ Could not compute ", which_end,
                      "adaper using heaviest path  method",
                      file=sys.stderr)
                print("\t/!\\ The resulting graph probably contains a loop.",
                      file=sys.stderr)
                adapters["heavy"][which_end] = ""

            # Exporting graph, if required
            if(args.export_graph is not None):
                if(v >= 1):
                    print("\tExporting assembly graph", file=print_dest)

                base = "_adapter_graph"
                if(".graphml" in args.export_graph):
                    base = "_" + os.path.basename(args.export_graph)
                else:
                    path = os.path.join(os.path.dirname(
                        args.export_graph), which_end + base + ".graphml")

                nx.write_graphml(g, path)
        else:
            unable_to_build = True

    if(not unable_to_build):
        # If we just need to print the adapter
        if(just_print):
            print_result(adapters, v, print_dest)
        # If we need to use them, getting back to porechop objects
        else:
            if(v >= 1):
                print("Building adapter object", file=print_dest)

            # Adding adapter if either end was found, for both method
            if(adapters["heavy"]['start'] or adapters["heavy"]['end']):
                adp.append(
                    Adapter("abinitio_heavy_adapter",
                            start_sequence=('abinitio_heavy_Top',
                                            adapters["heavy"]['start']),
                            end_sequence=('abinitio_heavy_Bottom',
                                          adapters["heavy"]['end'])))
            else:
                print("\t/!\\Heavy adapter was not added to adapter list",
                      file=sys.stderr)
            if(adapters["greedy"]['start'] or adapters["greedy"]['end']):
                adp.append(
                    Adapter("abinitio_greedy_adapter",
                            start_sequence=('abinitio_greedy_Top',
                                            adapters["greedy"]['start']),
                            end_sequence=('abinitio_greedy_Bottom',
                                          adapters["greedy"]['end'])))
            else:
                print("\t/!\\Greedy adapter was not added to adapter list",
                      file=sys.stderr)

            if(v >= 1):
                print("The inference of adapters sequence is done.",
                      file=print_dest)
    else:
        print("#################################", file=sys.stderr)
        print("WARNING - Unable to build graph  ", file=sys.stderr)
        print("Count file format may be invalid ", file=sys.stderr)
        print("   No adapter will be returned   ", file=sys.stderr)
        print("#################################", file=sys.stderr)
    return(adp)


##############################################################################
#                                  MAIN                                      #
##############################################################################


def execFindAdapt(args):
    """ Launch adaptFinder to perform the approximate kmers count and try
    to rebuild the adapter using different methods.
        The count will be stored as a temporary file in a ./tmp folder.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    # getting args
    fasta_file=args.input
    print_dest=args.print_dest
    v=args.verbosity

    # Temporary files prefix
    filename_pref="./tmp/temp_file_"
    tmpDir=os.path.dirname(filename_pref)

    # Creating tmp folder if not existing
    if(not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)

    # Searchng path to adaptFinder
    adapt_path=os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "adaptFinder")

    # using custom config file if specified
    conf_path=os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "ab_initio.config")
    conf_path=args.ab_initio_config if args.ab_initio_config else conf_path
    out_file_name=filename_pref + "approx_kmer_count"
    nb_thread=str(min(4, cpu_count()))
    if(v > 0):
        print("Using config file:" + conf_path, file = print_dest)

    # building command
    command=adapt_path + " " + fasta_file + " --config " + \
        conf_path + " -v " + str(v) + " -o " + out_file_name
    if(v > 0):
        print(command)

    # executing adapter finder algorithm
    try:
        subprocess.check_call(command.split())
    except SystemError as e:
        print(e, file = sys.stderr)
        sys.exit("ERROR: approximate k-mer count failed")

    # building the adapter from the count files generated
    adp=build_adapter(out_file_name, args)

    return(adp)


if __name__ == '__main__':
    adapt=execFindAdapt(sys.argv[1])
    for ad in adapt:
        print(ad.__dict__)
