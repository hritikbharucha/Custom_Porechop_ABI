"""
Ab-Initio wrapper created by Quentin Bonenfant
quentin.bonenfant@univ-lille.fr
quentin.bonenfant@gmail.com

https://github.com/qbonenfant
https://github.com/bonsai-team

This script evolved from the wrapper used to call the adaptFinder 
algorithm, which goal is to find adapter sequence k-mers.
We are able to infer adapter sequence from the reads instead of using
the adapter.py static database (unsupported since 2018) of native Porechop.
In order to change default parameters for adaptFinder, change the config
file porchop/abinitio.conf or pass a custom config file based on the same
syntax using -abc / --ab_initio_config

This works has been funded by the french National Ressearch Agency (ANR).

This file is part of Porechop_ABI. Porechop_ABI is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop_ABI. If not, see <http://www.gnu.org/licenses/>.
"""

import os
import subprocess
from multiprocessing import cpu_count
from .adapters import Adapter
from .consensus import find_consensus, build_consensus_graph
from collections import defaultdict as dd
from statistics import median, mean
import networkx as nx
import sys
from .arg_parser import get_arguments


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

    for meth in METHODS:
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

    out = print_dest
    for meth in METHODS:
        srt, end = "", ""
        if(meth in adapters.keys()):
            msg = meth
            if("start") in adapters[meth].keys():
                for k, v in adapters[meth]["start"].items():
                    srt += f"{k}: {v}\n"
            if("end") in adapters[meth].keys():
                for k, v in adapters[meth]["end"].items():
                    end += f"{k}: {v}\n"

            meth += " assembly method"
            meth = meth.capitalize()
            srt = "Start:\n" + srt
            end = "End:\n" + end

            print(msg, file=out)
            print(srt, file=out)
            print(end, file=out)


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


def execFindAdapt(args, out_file_name, mr, v, print_dest):

    #################################################################
    # PREPARING FILE SYSTEM
    #
    # Searching path to adaptFinder
    adapt_path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "adaptFinder")

    # Using custom config file if specified
    conf_path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "ab_initio.config")
    conf_path = args.ab_initio_config if args.ab_initio_config else conf_path

    # Input filename
    fasta_file = args.input

    # TODO: Make number of thread adjustable.
    nb_thread = str(min(4, cpu_count()))

    # Building command line for adaptFinder
    command = adapt_path + " " + fasta_file + \
        f" -v {v}" + \
        f" --config {conf_path}" + \
        f" -o {out_file_name}" + \
        f" -nt {nb_thread}"

    # if multi run, add mr flag
    if(mr > 1):
        command += f" -mr {mr}"

    if(v > 0):
        print("Using config file:" + conf_path, file=print_dest)
        print("Command line:\n", command,
              file=print_dest)

    try:
        # DEBUG
        # print("SUBPROCESS TRY", file=print_dest)
        subprocess.check_call(command.split())
        # subprocess.check_call(command.split(), stdout=print_dest)
        # subprocess.run(command.split(), check=True)
    except SystemError as e:
        print("\n#####################################",
              file=sys.stderr)
        print("ERROR: approximate k-mer count failed\n" + e,
              file=sys.stderr)
        print("#####################################\n",
              file=sys.stderr)
        sys.exit(1)


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


def make_adapter_object(adapters, v, print_dest):
    adp = []
    #################################################################
    # RESULT EXPORT
    #
    # If we just need to print the adapter
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
    return(adp)


def build_adapter(args, out_file_name, v, print_dest):
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

    # Launchiong adaptFinder in simple mode
    execFindAdapt(args, out_file_name, 1, v, print_dest)

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


def find_general_consensus(adpDict):
    """ Find a consensus for each method and end of adapter sequence.
    """


def consensus_adapter(args, prefix, v, print_dest):
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

    # if we only need to print the inferred adapter
    just_print = args.guess_adapter_only

    # Defining output filename basename
    out_file_name = prefix + "approx_kmer_count"

    #################################################################
    # CONSENSUS BUILDING
    #
    # Starting number of runs
    nb_run = args.multi_run
    # While we have no consensus, build new runs.
    consensus_found = False
    # First batch need to have 100% consensus, else we rerun
    first_batch_done = False
    # Unless an error occured.
    build_error = False

    while not consensus_found and not build_error:
        if(v > 0):
            if(not first_batch_done):
                print(f"Starting with a {nb_run} run batch.", file=print_dest)
            else:
                print(f"Following with a {nb_run} run batch.", file=print_dest)    

        # Launching adaptFinder in multi run mode
        execFindAdapt(args, out_file_name, nb_run, v, print_dest)

        for i in range(nb_run):
            current_adapter = build_adapter(args,
                                            out_file_name + f"_{i}",
                                            v,
                                            print_dest)
            # storing adapters
            insert_adapter_in_adpDict(current_adapter, adpDict)

        # If this is the first run:
        if(not first_batch_done):
            first_batch_done = True
            consensus_found = True
            # Testing if we have a consensus
            for m in adpDict.keys():
                adapters[m] = {}
                tmp = {}
                for e in adpDict[m].keys():
                    found = adpDict[m][e]
                    for a in found.keys():
                        # we only consider the consensus is found
                        # if all found adapter are the same in all
                        # runs
                        consensus_found &= found[a] == nb_run
                    best_adp = max(found.keys(), key=lambda x: found[x])
                    tmp[e] = best_adp

                adapters[m] = tmp

        # If we are at the second run, we need to find a consensus
        else:
            consensus_found = True
            # for each method and ends
            for m in adpDict.keys():
                adapters[m] = {}
                tmp = {}
                for e in adpDict[m].keys():
                    g = build_consensus_graph(adpDict[m][e])
                    path = find_consensus(g, nb_run)
                    tmp[e] = concat_path(path)
                adapters[m] = tmp

        # Do we need to rerun ?
        if(not consensus_found):
            # If a 100% consensus is not found at first,
            # we increase the number of run and add those results
            # to current adapter collection (default: 20)
            nb_run = args.consensus_run
            out_file_name += "_sup"
            print("/!\\\tMore runs are required to build consensus.",
                  file=sys.stderr)
            print(f"/!\\\tSetting number of run to {nb_run}.",
                  file=sys.stderr)
            print("/!\\\tCurrent adapter distribution:",
                  file=sys.stderr)
            print_adapter_dict(adpDict, sys.stderr)

    if(consensus_found):
        if(v > 0):
            print("consensus step done", file=print_dest)
        if(args.export_consensus):
            try:
                out = open(args.export_consensus, "at")
            except FileNotFoundError:
                print("Could not export consensus file to", file=sys.stderr)
                print(args.export_consensus)
            else:
                print_adapter_dict(adpDict, out)
                out.close()

    #################################################################
    # RESULT EXPORT
    #
    # If we just need to print the adapter
    if(just_print):
        print_result(adapters, v, print_dest)

    # If we need to use them, we build porechop Adapter objects
    else:
        adp = make_adapter_object(adapters, v, print_dest)
    return adp


##############################################################################
#                                  MAIN                                      #
##############################################################################


def launch_ab_initio(args):
    """ Launch adaptFinder to perform the approximate kmers count and try
    to rebuild the adapter using different methods.
        The count will be stored as a temporary file in a ./tmp folder.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    # getting args
    print_dest = args.print_dest
    v = args.verbosity

    # Temporary files prefix
    filename_pref = "./tmp/temp_file_"
    tmpDir = os.path.dirname(filename_pref)

    # Creating tmp folder if not existing
    if(not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)

    # If multiple run are required:
    if(args.multi_run > 1):
        # executing consensus adapter inference algorithm
        adp = consensus_adapter(args, filename_pref, v, print_dest)
    else:
        # Or regular adapterinference algorithm
        out_file_name = filename_pref + "approx_kmer_count"
        adapters = build_adapter(args, out_file_name, v, print_dest)
        adp = make_adapter_object(adapters, v, print_dest)
    return(adp)


if __name__ == '__main__':
    args = get_arguments
    adapt = launch_ab_initio(args)
    for ad in adapt:
        print(ad.__dict__)
