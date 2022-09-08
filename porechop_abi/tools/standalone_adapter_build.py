"""
Standalone script for adapter reconstruction from k-mer count file.
This script is more is able to rebuild consensus run using a folder
of count file (if you kept them). It is an heavily modified version
of the script that ships with adpatFinder, and is tailored for
Porechop_ABI count file.

"""

import os
import argparse
from consensus import find_consensus, build_consensus_graph
from collections import defaultdict as dd
from statistics import median, mean
import shutil
import networkx as nx
import sys


##############################################################################
#                                 CONSTANTS                                  #
##############################################################################

CUT_RATIO = 0.05  # How much random variation is allowed for drop cut.
METHODS = ["greedy", "heavy"]  # Reconstruction methods.
ENDS = ["start", "end"]        # Possible ends key.


##############################################################################
#                              UTILITY FUNCTIONS                             #
##############################################################################


def haveOverlap(seq1, seq2):
    """Check if the sequence 1 is a prefix of sequence 2
    @param first sequence
    @param second sequence
    @return seq1 and seq2 merged in one sequence on the overlap.
    """
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


def get_weight(g, path):
    """Compute the weight of a path in an Nx graph
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
    s_cut = len(start_adp_data) - 1
    if(smz):
        s_cut = min(max(smz), s_cut)

    # Returning the cut path
    return(start_adp_data[0:s_cut])


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
    e_cut = 0
    if(emz):
        e_cut = max(0, min(emz) + 1)

    # Returning the cut path
    return(end_adp_data[e_cut:])


def sliding_median(counts, w=7):
    """ Compute the median of a of each count value using
    a sliding window approach.
    @param counts a list of integer
    @param w the window size
    @return the median values, in a list.
    """
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
    """ Compute the mean of a of each count value using
    a sliding window approach.
    @param counts a list of integer
    @param w the window size
    @return the mean values, in a list.
    """
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
    """ Compute the difference between consecutives values in a list
    @param counts a list of integer
    @return The difference between a value and it's neighbourg.
    """
    lc = len(counts)
    results = [0]
    for i in range(1, lc, 1):
        val = counts[i] - counts[i - 1]
        results.append(val)
    return(results)


def find_median_zone(counts, epsilon):
    """ Find the position in the counts list where the
    counts raise above a threshold.
    @param counts the count list
    @param epsilon margin of acceptable error
    """
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


def print_result(adapters, v, print_dest=sys.stdout):
    """Display the result, adjusting format depending on selected verbosity
    @param an adapter dictionnary, using appropriate methods as keys
    @param v, the verbosity level
    @param print destination, can be custom, but stdout is default
    """

    global METHODS

    out = print_dest

    if(v > 0):
        print("\n\nINFERRED ADAPTERS:\n",
              file=out)
    else:
        out = sys.stdout

    for meth in sorted(METHODS, reverse=True):
        msg = meth
        srt = adapters[meth]["start"]
        end = adapters[meth]["end"]
        if(v >= 1):
            meth += " assembly method"
            meth = meth.capitalize()
            srt = "Start:\t" + srt
            end = "End:\t" + end
        print(msg, file=out)
        print(srt, file=out)
        print(end, file=out)


def print_adapter_dict(adapters, print_dest=sys.stderr):
    """Display the content of adapter counting dictionnaries.
    @param an adapter dictionnary, using appropriate methods as keys
    @param print destination, can be custom, but stdout is default
    """

    global METHODS

    for meth in sorted(METHODS, reverse=True):
        srt, end = "", ""
        if(meth in adapters.keys()):
            msg = meth
            if("start") in adapters[meth].keys():
                for k, v in adapters[meth]["start"].items():
                    srt += f"{k}: {v}\n"
            if("end") in adapters[meth].keys():
                for k, v in adapters[meth]["end"].items():
                    end += f"{k}: {v}\n"

            msg += " assembly method"
            msg = msg.capitalize()
            srt = "Start:\n" + srt
            end = "End:\n" + end

            print(msg, file=print_dest)
            print(srt, file=print_dest)
            print(end, file=print_dest)


def clean_input_filename(infile):
    if(".start" == infile[-6:] or ".end" == infile[-4:]):
            infile = ".".join(infile.split(".")[:-1])
    return(infile)




##############################################################################
#                                 ASSEMBLY                                   #
##############################################################################

def greedy_assembl(g):
    """Greedy assembly method to compute the adapter using the graph as input.
    @param the De Bruijn graph of kmers
    @return the longest debruijn sequence starting by the most frequent kmer
    """
    start = max(g.nodes, key=lambda x: g.nodes[x]["weight"])
    path = [start]

    right_node = start
    left_node = start
    while(left_node or right_node):

        # forward extension
        if(right_node):
            r_list = [el for el in g.successors(right_node) if el not in path]
            if(not r_list):
                right_node = None
            else:
                right_node = max(r_list, key=lambda x: g.nodes[x]["weight"])
                path.append(right_node)

        # reverse extension
        if(left_node):
            l_list = [el for el in g.predecessors(left_node) if el not in path]

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


def generic_build(method_fct, name, g, adapters, which_end, v, w, print_dest):
    """Generic function used to build adapter from a path in the graph g.
    This function need a path building function (which returns a list of
    consecutive overlapping k-mers). The resulting adapter is then
    added in the supplied adapters dictionnary for further processing.
    """
    if(v >= 1):
        print(f"\tBuilding {name} {which_end} adapter", file=print_dest)

    try:
        method_p = method_fct(g)
    except(nx.exception.NetworkXUnfeasible,
           nx.exception.NetworkXNotImplemented,
           nx.exception.HasACycle) as nxE:

        print(f"\t/!\\Could not compute {which_end} adaper with {name} method",
              file=sys.stderr)
        print("\t/!\\ The resulting graph probably contains a loop.",
              file=sys.stderr)
        print("\t/!\\ For more details, read the error report below;\n",
              file=sys.stderr)
        print(nxE, file=sys.stderr)
        adapters[name][which_end] = ""

    else:
        cut_method_p = []
        if(which_end == "start"):
            cut_method_p = start_cut(method_p, g, w)
        elif(which_end == "end"):
            cut_method_p = end_cut(method_p, g, w)
        # If the method path is not empty, concat it in a string
        if(cut_method_p):
            adapters[name][which_end] = concat_path(cut_method_p)
        else:
            adapters[name][which_end] = ""



def build_adapter(args, out_file_name, v, print_dest, export_graph=""):
    """Building adapters from counts using different method:
        - Building a Debruijn graph and searching heaviest path
        - Greedy assembly based on kmer rank (most frequent first)
    Porechop can still properly detect adapters if there is a slight
    difference  since it uses a quite low identity threshold of 75%.
    @param A De Bruijn graph with kmer count as weight
    @param The argparse argument dictionnary from porechop.py
    @return A list of Adapter objects containing infered adapters.
    """
    w = args.window_size
    adapters = dd(dict)  # temporary adapter string dict

    # Failsafe if the graphs for the adapter can not be built
    unable_to_build = False

    for which_end in ENDS:
        if(v >= 1):
            print("Assembling " + which_end + " adapters", file=print_dest)

        # Building graph
        g = build_graph(out_file_name + "." + which_end)

        # Checking graph is properly build
        if(not nx.is_empty(g)):

            # preping for automatic anotation
            nx.set_node_attributes(g, "", "path")

            # Greedy Adapter
            generic_build(greedy_path, "greedy",
                          g,
                          adapters,
                          which_end,
                          v,
                          w,
                          print_dest)
            # Heavy adapter
            generic_build(heavy_path, "heavy",
                          g,
                          adapters,
                          which_end,
                          v,
                          w,
                          print_dest)

            # Exporting graph, if required
            if(export_graph):
                if(v >= 1):
                    print("\tExporting assembly graph", file=print_dest)
                base = os.path.basename(export_graph)
                base = base[:-8] if ".graphml" == base[-8:] else base
                path = os.path.dirname(export_graph)
                g_file = os.path.join(path, f"{base}_{which_end}.graphml")
                nx.write_graphml(g, g_file)
        else:
            unable_to_build = True

    if(unable_to_build):
        print("#################################", file=sys.stderr)
        print("WARNING - Unable to build graph  ", file=sys.stderr)
        print("Count file format may be invalid ", file=sys.stderr)
        print("   No adapter will be returned   ", file=sys.stderr)
        print("#################################", file=sys.stderr)

    return adapters

##############################################################################
#                            Multi Run / Consensus                           #
##############################################################################


def insert_adapter_in_adpDict(adapter, adpDict):
    """Insert an adapter in an adapter collection dictionnary
    @param adapter An nested adapter dictionnary (method used -> start/end)
    @param adpDict A collection of previous adapter (includes adp frequency)
    """
    for method in adapter.keys():
        adp_ends = adapter[method]
        for which_end in adp_ends.keys():
            if(method not in adpDict.keys()):
                adpDict[method] = {}
            if(which_end not in adpDict[method].keys()):
                adpDict[method][which_end] = {}
            the_adapter = adapter[method][which_end]
            if(the_adapter in adpDict[method][which_end].keys()):
                adpDict[method][which_end][the_adapter] += 1
            else:
                adpDict[method][which_end][the_adapter] = 1



def rebuild_consensus_adapter(args, folder, v, print_dest):
    """ Try to find the fittest adapter for the dataset by running multiple
    adapter reconstruction and building a consensus.
    We first try 10 runs, and keep the adapter if a perfectly stable consensus
    is found immediatly. Otherwise, 20 other runs are launched and the adapter
    sequence is build from a consensus of the 30 runs.
    """

    # Final adapters placeholders
    adp = []
    adapters = {}

    # Consensus adapter placeholder
    adpDict = {}

    # Defining output filename for graphs
    g_directory = os.path.dirname(args.export_graph) if(args.export_graph) else ""

    # If we have to keep initial / supplementary run separated:
    # build two structure to store them separately
    if(args.split_run_types):
        initial_adp = {}
        supplem_adp = {}
    

    #################################################################
    # CONSENSUS BUILDING
    #

    file_set = set(clean_input_filename(el) for el in  os.listdir(folder))
    nb_run = len(file_set)
    #while not consensus_found and not build_error:
    for prefix in list(file_set):
        eg = ""
        if(g_directory):
            eg = os.path.join(g_directory, prefix)

        out_file_name = os.path.join(folder, prefix)
        current_adapter = build_adapter(args,
                                        out_file_name,
                                        v,
                                        print_dest,
                                        eg)
        # storing adapters
        insert_adapter_in_adpDict(current_adapter, adpDict)
        # storing separated too if wanted
        if(args.split_run_types):
            if("_sup_" in out_file_name):
                insert_adapter_in_adpDict(current_adapter, supplem_adp)
            else:
                insert_adapter_in_adpDict(current_adapter, initial_adp)

    # Building the consensus for each method and ends
    for m in adpDict.keys():
        adapters[m] = {}
        tmp = {}
        for e in adpDict[m].keys():
            g = build_consensus_graph(adpDict[m][e])
            path = find_consensus(g, nb_run)
            tmp[e] = concat_path(path)
        adapters[m] = tmp
    
    
    #################################################################
    # RESULT EXPORT
    #

    out = open("consensus_export.txt", "w")
    # Print all adapter found to a file
    if(args.split_run_types):
        print("##############################", file=out)
        print("Initial runs", file=out)
        print_adapter_dict(initial_adp, out)
        print("##############################", file=out)
        print("Supplementary runs", file=out)
        print_adapter_dict(supplem_adp, out)
    else:
        print_adapter_dict(adpDict, out)
    out.close()

    # Print resulting adapter
    print_result(adapters, v, print_dest)



##############################################################################
#                            ARGUMENT PARSING                                #
##############################################################################
class MyHelpFormatter(argparse.HelpFormatter):

    """
    This is a custom formatter class for argparse. It allows for some custom
    formatting, in particular for the help texts with multiple options
    (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if action.default != argparse.SUPPRESS and \
           'default' not in help_text.lower() and \
           action.default is not None:

            help_text += ' (default: ' + str(action.default) + ')'
        return help_text


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='Build adapter'
                                     'A tool for building adapters in Oxford '
                                     'Nanopore reads from k-mer count.',
                                     formatter_class=MyHelpFormatter,
                                     add_help=False)

    # Either a single file or a folder
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', default = "",
                            help='Path to k-mer count file or file name prefix\
                            (without the ".start" or ".end)".')
    input_group.add_argument('-mrf', '--multirun_folder', type=str, default="",
                            help='Path to the count file folder if you want a'
                                 'consensus over multiple runs')

    trim_group = parser.add_argument_group('Trimming')
    trim_group.add_argument('-w', '--window_size', type=int, default=3,
                           help='Size of the sliding window for smoothing')


    opt_group = parser.add_argument_group('Options')
    opt_group.add_argument('-v', '--verbosity', type=int, default=0,
                            help='Level of info display 0 = minimal, 1 = std,'
                                 '2 = a lot.')
    
    opt_group.add_argument('-e', '--export_graph', type=str,
                           help='Path to export the graph used for assembly\
                           (.graphml format), if you want to keep it')

    opt_group.add_argument('-srt', '--split_run_types', action='store_true',
                           default=False,
                           help='Multi run count folder only. When exporting\
                           the consensus sequences, split the initial run \
                           files and the supplementary runs when exporting\
                           adapters before consensus.')


    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help',
                           default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    

    args = parser.parse_args()
    return args



##############################################################################
#                                  MAIN                                      #
##############################################################################


def main(args):
    """ Use the supplied approximate kmers count to try to rebuild the adapter
        using different methods.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    # getting args
    print_dest = sys.stdout
    v = args.verbosity

    # If multiple mode are required:
    if(args.multirun_folder):
        # executing consensus adapter inference algorithm
        folder = args.multirun_folder
        rebuild_consensus_adapter(args, folder, v, print_dest)
    else:
        # In single file mode, we reshape the count file name to only keep the prefix.
        infile = clean_input_filename(args.input)
        # Then we use regular adapter building method
        adapters = build_adapter(args, infile, v, print_dest, args.export_graph)
        print_result(adapters, v, print_dest)



if __name__ == '__main__':
    args = get_arguments()
    main(args)
