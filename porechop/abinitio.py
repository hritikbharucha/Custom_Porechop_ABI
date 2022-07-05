"""
Ab-Initio inference of ONT adapter main file

https://github.com/qbonenfant
https://github.com/bonsai-team

This script evolved from the wrapper used to call the approx_counter 
algorithm, which end goal is to find adapter sequence k-mers.
We are able to infer adapter sequence from the reads instead of using
the adapter.py static database (unsupported since 2018) of native Porechop.
In order to change default parameters for approx_counter, change the config
file porchop/abinitio.conf or pass a custom config file based on the same
syntax using -abc / --ab_initio_config

This works has been funded by the french National Ressearch Agency (ANR).

This file is part of Porechop_ABI. Porechop_ABI is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop_ABI is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop_ABI. If not, see <http://www.gnu.org/licenses/>.
Author: Quentin Bonenfant (quentin.bonenfant@gmail.com)
"""

import os
from multiprocessing import cpu_count
from .adapters import Adapter
from collections import defaultdict as dd
import networkx as nx
import sys
from .arg_parser import get_arguments
from .drop_cut import *
from .consensus import *

##############################################################################
#                                 CONSTANTS                                  #
##############################################################################

METHODS = ["greedy", "heavy"]  # Reconstruction methods.
ENDS = ["start", "end"]        # Possible ends key.
LOW_FREQ_THRESHOLD = 10        # K-mer low frequency treshold.


##############################################################################
#                              UTILITY FUNCTIONS                             #
##############################################################################


def warning(msg):
    """
    Print a message to stderr in a formatted way
    @param msg: the message to print as warning (str)
    """
    print(f"/!\\\t {msg}", file=sys.stderr)


def low_freq_warning(end, amount):
    
    warning(f"The most frequent kmer has been found in less"+
            f" than {LOW_FREQ_THRESHOLD}% of the reads {end}s"+
            f" after approximate count ({amount}%)")
    warning("It could mean this file is already trimmed or the sample do not"+
             " contains detectable adapters.")


def print_consensus_result(consensus_adapters, v, print_dest=sys.stderr):
    """Display the result from consensus runs. The format is different 
    beacause methods are pooled in consensus runs, so a method oriented 
    display no longer makes sens.
    @param an adapter dictionnary, using appropriate ends as keys
    @param v, the verbosity level
    @param print destination, can be custom, but stdout is default
    """

    if(v > 0):
        print("\n\nCONSENSUS ADAPTERS:\n",
              file=print_dest)

    for end in ENDS:
        print(end.capitalize(), file=print_dest)
        sorted_consensus = sorted(consensus_adapters[end].keys(),
                                  key=lambda x: int(x.split("_")[1]))
        for name in sorted_consensus:
            adp = consensus_adapters[end][name]
            print(name, file=print_dest)
            print(adp, file=print_dest)


def print_result(adapters, v, print_dest=sys.stderr):
    """Display the result in the legacy format (method->ends)
    using the updated data structure (end -> methods)
    @param an adapter dictionnary, using appropriate methods as keys
    @param v, the verbosity level
    @param print destination, can be custom, but stdout is default
    """
    out = print_dest

    if(v > 0):
        print("\n\nINFERRED ADAPTERS:\n",
              file=print_dest)

    # Getting methods from the dict directly.
    methods = sorted(adapters[ENDS[0]].keys())

    for meth in methods:
        msg = meth
        srt = adapters["start"][meth]
        end = adapters["end"][meth]
        if(v > 0):
            meth = meth.capitalize()
            srt = "Start:\t" + srt
            end = "End:\t" + end
        print(msg, file=print_dest)
        print(srt, file=print_dest)
        print(end, file=print_dest)


def print_adapter_dict(adapters, print_dest=sys.stderr):
    """Display the content of adapter counting dictionnaries.
    using the updated dict format.
    @param an adapter dictionnary, using appropriate methods as keys
    @param print destination, can be custom, but stderr is default
    """

    methods = sorted(adapters[ENDS[0]].keys())

    out = print_dest
    for meth in methods:
        srt, end = "", ""
        msg = meth
        for k, v in adapters["start"][meth].items():
            srt += f"{k}: {v}\n"
        for k, v in adapters["end"][meth].items():
            end += f"{k}: {v}\n"

        meth += " assembly method"
        meth = meth.capitalize()
        srt = "Start:\n" + srt
        end = "End:\n" + end

        print(msg, file=print_dest)
        print(srt, file=print_dest)
        print(end, file=print_dest)

def parse_conf(conf_file):
    """
    Parse the configuration file passed to approxCounter.
    It allows Porechop_ABI to know which parameters are 
    used during the approximate k-mer count.
    @param The path to the configuration file (str)
    @return A dictionnary containing the parameters.
    """

    # Setting some default values:
    params = {'k' : 16,      # kmer size, 2<= k <= 32
              'sl' : 100,     # sequence sampling size
              'sn' : 40000,   # number of sequence sampled
              'limit' : 500,  # number of kmers to keep.
              'lc' :    1.2   # low complexity filter threshold.
    }
    try:
        assert(os.path.isfile(conf_file))
    except AssertionError:
        warning("Invalide configuration file.")
        warning(conf_file)
        sys.exit(1)
    else:
        with open(conf_file) as f:
            for line in f:
                line = line.rstrip("\n").replace(" ","")
                data = line.split("=")
                if(line and line[0]!="#" and len(data) == 2):
                    p, v = data
                    # int values
                    if(v.isdigit()):
                        params[p] = int(v)
                    # float values
                    elif(v.replace(".","",1).isdigit()):
                        params[p] = float(v)
                    # string values
                    else:
                        params[p] = v
    return(params)




##############################################################################
#                           ADAPTER SEQUENCE ASSEMBLY                        #
##############################################################################

def greedy_assembly(g):
    """Greedy assembly method to compute the adapter using the graph as input.
    @param  g: A NetworkX DiGraph (graph)
    @return path: The longest extendable sequence from the most frequent kmer.
    """

    # Finding the best node in the graph.
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

    @param G: A NetworkX DiGraph (graph)
    @return path: the heaviest path in the graph (list)

    Comment:    
    This is a modified version of the dag_longest_path from Networkx
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
    """ Searching the heaviest path between all source and target nodes.
        Even if longer path tend to be heavier, very heavy short path can
        also be selected.
        @param g: A NetworkX DiGraph (graph)
        @return hv_path : The heaviest path in the graph using node weight (list)
    """

    hv_path = dag_heaviest_path(g)

    # Annotating graph
    for n in hv_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",heavy").lstrip(",")

    return(hv_path)


def greedy_path(g):
    """Greedy assembly of the adapter from the graph,
    starting from the most frequent
    @param  g: A NetworkX DiGraph (graph)
    @return path: The longest extendable sequence from the most frequent kmer.
    """

    gd_path = greedy_assembly(g)
    # annotating graph
    for n in gd_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",greedy").lstrip(",")
    return(gd_path)


##############################################################################
#                              GRAPH FUNCTIONS                               #
##############################################################################

def haveOverlap(seq1, seq2):
    """Check if the sequence 1 is a k-1 prefix of sequence 2
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
    @param g: The Networkx DiDraph (graph)
    @param path: The path (list)
    @return total: The total weight of this path (int)
    """
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)


def concat_path(path):
    """Concat the kmers of a path into a single sequence
    @param The path as a k-mer list (list)
    @return the full sequence. (string)
    """
    return(path[0][:-1] + "".join(el[-1] for el in path))



def build_graph(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
       @param count_file : the path to the file to build the graph (string)
       @return g: A Networkx DiGraph containing all weighted k-mers (graph)
    """
    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # Try to build the graph.
    try:
        with open(count_file, 'r') as f:
            for line in f:
                km, nb = line.rstrip("\n").split("\t")
                kmer_count[km] = int(nb)
                kmer_list.append(km)
                g.add_node(km, weight=int(nb))  # , path = "" )

    # If the file can not be opened.
    except FileNotFoundError:
        warning("Unable to open k-mer count file:")
        warning(count_file)
        warning("Either the file was moved, deleted, or filename is invalid.")
    
    # in case we try to convert to int something that is not a number
    except TypeError:
        warning("Something is wrong with the k-mer count file format:")
        warning(count_file)
        warning("Each line need to be in this format: <k-mer> \\t <count>")

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


def execApproxCounter(args, out_file_name, conf_file, mr, v, print_dest=sys.stderr):
    """ Prepare command and execute the approximate k-mer counting 
    program using subprocess.
    @param args: The arguments passed to Porechop (argparse arguments)
    @param out_file_name : The prefix for the k-mer count output file.
    @param mr : number of run to perform (for multi run). default is 1 (int)
    @param v: verbosity level. Default is 1. (int)
    @param print_dest : File to print to, std out by default (file stream)

    Comment: This function returns nothing, it either generate a k-mer counting file
    or stops the whole program.
    """

    # Searching path to approx_counter
    prog_path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "approx_counter")

    # Input filename
    fasta_file = args.input

    # TODO: Manage number of threads in a better way.
    nb_thread = str(min(args.threads, cpu_count()))

    # Building command line for approx_counter
    command = prog_path + " " + fasta_file + \
        f" -v {v}" + \
        f" --config {conf_file}" + \
        f" -o {out_file_name}" + \
        f" -nt {nb_thread}"

    # If multi run, add mr flag to command line
    if(mr > 1):
        command += f" -mr {mr}"

    if(v > 0):
        print("Using config file:" + conf_file, file=print_dest)
        print("Command line:\n", command,
              file=print_dest)

    # Finally, Running subprocess to count the k-mers
    try:
        subprocess.check_call(command.split())
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
    @param method_fct: Handle to the function used to build the path. (python function)
    @param name : The name of the method used (string)
    @param g :  A NetworkX DiGraph (graph)
    @param adapters: An adapter dictionnary to store results. (defaultdict)
    @param which_end: Name of the sequence end currently processed: start or end (string)
    @param v : verbosity level. Default is 1. (int)
    @param w : Window size for the drop cut mean smoothing
    @param print_dest : File to print to, std out by default (file stream)
    """
    if(v >= 1):
        print(f"\tBuilding {name} {which_end} adapter", file=print_dest)

    try:
        method_p = method_fct(g)
    except(nx.exception.NetworkXUnfeasible,
           nx.exception.NetworkXNotImplemented,
           nx.exception.HasACycle) as nxE:

        warning(f"Could not compute {which_end} adaper with {name} method")
        warning("The resulting graph probably contains a loop.")
        warning("For more details, read the error report below:")
        warning(nxE)
        warning("\n")
        adapters[which_end][name] = ""

    else:
        cut_method_p = []
        if(which_end == "start"):
            cut_method_p = start_cut(method_p, g, w)
        elif(which_end == "end"):
            cut_method_p = end_cut(method_p, g, w)
        # If the method path is not empty, concat it in a string
        if(cut_method_p):
            adapters[which_end][name] = concat_path(cut_method_p)
        else:
            warning(f"Something went wrong for {which_end} \
                  adaper with {name} method")
            warning("The returned graph path is empty")
            adapters[which_end][name] = ""


def make_adapter_objects(adapters, v, print_dest):
    adp = []
    # If we just need to print the adapter
    if(v >= 1):
        print("Building adapter object", file=print_dest)

    methods = sorted(adapters[ENDS[0]].keys())

    for meth in methods:
        adp.append(
            Adapter(f"abinitio_{meth}_adapter",
                    start_sequence=(f'abinitio_{meth}_Top',
                                    adapters['start'][meth]),
                    end_sequence=(f'abinitio_{meth}_Bottom',
                                  adapters['end'][meth])))

    if(v >= 1):
        print("The inference of adapters sequence is done.",
              file=print_dest)
    return(adp)


def build_adapter(args, out_file_name, sn, v, print_dest):
    """Building adapters from counts using different method:
        - Building a Debruijn graph and searching heaviest path
        - Greedy assembly based on kmer rank (most frequent first)
    Porechop can still properly detect adapters if there is a slight
    difference  since it uses a quite low identity threshold of 75%.
    @param A De Bruijn graph with kmer count as weight
    @param The argparse argument dictionnary from porechop.py
    @return A list of Adapter objects containing infered adapters.
    """

    # window size for the count smoothing at the drop cut stage.
    w = args.window_size

    # Keeping track of the best k-mer frequency
    best_freq = []

    # Placeholder adapter string dict
    adapters = dd(dict)  

    # Failsafe if the graphs for the adapter can not be built
    unable_to_build = False

    for which_end in ENDS:
        if(v >= 1):
            print("Assembling " + which_end + " adapters", file=print_dest)

        # Building graph
        g = build_graph(out_file_name + "." + which_end)

        # Checking graph is properly build
        if(not nx.is_empty(g)):

            # Computing best k-mer frequency to decide
            # if a warning is needed or not
            max_count = max(g.nodes[n]["weight"] for n in g.nodes)
            best_freq.append(max_count * 100 / sn)
            
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
            break

    if(unable_to_build):
        print("#################################", file=sys.stderr)
        print("ERROR  -  Unable to build graph  ", file=sys.stderr)
        print("Count file format may be invalid ", file=sys.stderr)
        print("#################################", file=sys.stderr)
        sys.exit(1)

    return (adapters, best_freq)

##############################################################################
#                            Multi Run / Consensus                           #
##############################################################################


def insert_adapter_in_adp_count(adapter, adp_count):
    """Insert an adapter in an adapter collection dictionnary
    @param adapter A nested adapter dictionnary (start/end -> method used)
    @param adp_count A collection of previous adapter (including frequency)
    """
    for which_end in adapter.keys():
        methods = adapter[which_end].keys()
        for method in methods:    
            # adding ends dictionnary to adp_count
            if(which_end not in adp_count.keys()):
                adp_count[which_end] = {}
            # adding method sub_dictionnary to adp_count
            if(method not in adp_count[which_end].keys()):
                adp_count[which_end][method] = dd(int)

            # updating adapter count
            the_adapter = adapter[which_end][method]
            adp_count[which_end][method][the_adapter] += 1


def make_consensus_adapter_objects(adapters, v, print_dest):
    """ Build Adapter object with only start or end sequence.
    Pairing consensus sequences do not make sense since they
    come from different samples.
    @param adapters: The nested dictionnary with adapter sequences
    @param v: verbosity level (int)
    @param print_dest: the file to print the info (stream) 
    """
    adp = []
    if(v >= 1):
        print("Building consensus adapter objects", file=print_dest)

    # start
    for consensus_name in adapters["start"].keys():
        adp.append(
            Adapter(f"{consensus_name}_adapter",
                    start_sequence=(f'{consensus_name}_Top',
                                    adapters['start'][consensus_name])))
    # end
    for consensus_name in adapters["end"].keys():
        adp.append(
            Adapter(f"{consensus_name}_adapter",
                    end_sequence=(f'{consensus_name}_Bottom',
                                  adapters['end'][consensus_name])))

    if(v >= 1):
        print("The inference of adapters sequence is done.",
              file=print_dest)
    return(adp)


def consensus_adapter(args, prefix, conf_file, v, print_dest):
    """ Try to find the fittest adapter for the dataset by running multiple
    adapter reconstruction and building a consensus.
    We first try 10 runs, and keep the adapter if a perfectly stable consensus
    is found immediatly. Otherwise, 20 other runs are launched and the adapter
    sequence is build from a consensus of the 30 runs.
    """

    # Final adapter object list
    adp = [] 
    # Adapters placeholders
    adapters = dd(dict)

    # Consensus adapter counter
    adp_count = {}

    # if we only need to print the inferred adapter
    just_print = args.guess_adapter_only

    # Defining output filename basename
    out_file_name = prefix + "approx_kmer_count"

    # getting params sent to approxCounter
    params = parse_conf(conf_file)

    # Using the sample size to determine if the k-mer frequency is small.
    sn = params["sn"]

    #################################################################
    # CONSENSUS BUILDING
    #
    # Starting number of runs
    nb_run = args.multi_run
    total_run = nb_run
    # While we have no consensus, build new runs.
    consensus_found = False
    # First batch need to have 100% consensus, else we rerun
    first_batch_done = False
    # Unless an error occured.
    build_error = False

    # looking for the best k-mer frequency for each runs
    run_best_km_freq = []
    while not consensus_found and not build_error:
        if(v > 0):
            if(not first_batch_done):
                print(f"Starting with a {nb_run} run batch.", file=print_dest)
            else:
                print(f"Following with a {nb_run} run batch.", file=print_dest)    

        # Launching approx_counter in multi run mode
        execApproxCounter(args, out_file_name, conf_file, nb_run, v, print_dest)

        # Parsing individual files
        for i in range(nb_run):
            current_adp, best_freq = build_adapter(args,
                                                   out_file_name + f"_{i}",
                                                   sn,
                                                   v,
                                                   print_dest)
            # storing and counting adapters
            insert_adapter_in_adp_count(current_adp, adp_count)
            # storing best k-mer frequency
            run_best_km_freq.append(best_freq)

        # pooling methods for consensus mode.
        pooled_name = "pooled"
        pooled = pool_methods(adp_count, pooled_name)
        # If this is the first run:
        if(not first_batch_done):
            first_batch_done = True
            consensus_found = True
            # Do we have a perfect consensus ?
            for end in ENDS:
                adp_list, counts = zip(*pooled[end][pooled_name].items())
                # if we have more than one adapter in the pooled dict, then no
                if(len(adp_list) != 1):
                    consensus_found = False
                else:
                    adapters[end]["Consensus_1_100%"] = adp_list[0]

        # If we are at the second run, we need to find a consensus
        else:
            consensus_found = True
            # building the consensus dictionnary from the pooled adapter dict.
            adapters = build_consensus_adapter_dict(args, pooled, total_run * 2)
            # checking adapter list is empty
            if(len(adapters) == 0):
                warning()


        # Do we need to rerun ?
        if(not consensus_found):
            # If a 100% consensus is not found at first,
            # we increase the number of run and add those results
            # to current adapter collection (default: 20)
            nb_run = args.consensus_run
            total_run += nb_run
            out_file_name += "_sup"
            warning("More runs are required to build consensus.")
            warning(f"Setting number of run to {nb_run}.")
            # warning("Current adapter distribution:")
            # print_adapter_dict(adp_count, print_dest)

        # In case we want to export temp adapters used for the consenus
        if(args.export_consensus):
            try:
                out = open(args.export_consensus, "at")
            except FileNotFoundError:
                print("Could not export consensus file to", file=sys.stderr)
                print(args.export_consensus)
            else:
                print_adapter_dict(adp_count, out)
                out.close()

    if(v > 0):
        print("consensus step done", file=print_dest)

    #################################################################
    # RESULT EXPORT
    #
    # Are their any frequency warnings ?
    for i in range(len(ENDS)):
        lowest = 100
        nb_below = 0
        for freq in run_best_km_freq:
            lowest = min(lowest, freq[i])
            nb_below += 1 if freq[i] < LOW_FREQ_THRESHOLD else 0
        if(nb_below):
            low_freq_warning(ENDS[i], lowest)
            warning(f"Warning triggered by {nb_below}/{total_run} runs")

    # If we just need to print the adapter
    if(just_print):
        print_consensus_result(adapters, v, sys.stdout)

    # If we need to use them, we build porechop Adapter objects
    else:
        adp = make_consensus_adapter_objects(adapters, v, print_dest)
    return adp


##############################################################################
#                                  MAIN                                      #
##############################################################################


def launch_ab_initio(args):
    """ Launch approx_counter to perform the approximate kmers count and try
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

    # temp adptater object placeholder
    adp = []
    # Creating tmp folder if not existing
    if(not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)

    # Using custom config file if specified
    conf_file = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "ab_initio.config")
    conf_file = args.ab_initio_config if args.ab_initio_config else conf_file

    # If multiple run are required:
    if(args.multi_run > 1):
        # executing consensus adapter inference algorithm
        adp = consensus_adapter(args, filename_pref, conf_file, v, print_dest)
    else:
        # Or regular adapter inference algorithm
        out_file_name = filename_pref + "approx_kmer_count"

        # Launching approx_counter in simple mode
        execApproxCounter(args, out_file_name, conf_file, 1, v, print_dest)

        # getting params sent to approxCounter
        params = parse_conf(conf_file)

        # Using the sample size to determine if the k-mer frequency is small.
        sn = params["sn"]

        # Finally building  the adapter
        adapters , best_freq = build_adapter(args, out_file_name, sn, v, print_dest)
        
        # Low frequency warnings
        for i in range(len(ENDS)):
            if(best_freq[i]< LOW_FREQ_THRESHOLD):
                low_freq_warning(ENDS[i], best_freq[i])

        # If we need to just print the adapter, do it
        if(args.guess_adapter_only):
            print_result(adapters, v, sys.stdout)
        # else export adapter as Adapter object
        else:
            adp = make_adapter_objects(adapters, v, print_dest)
    return(adp)


if __name__ == '__main__':
    args = get_arguments()
    adapt = launch_ab_initio(args)
    for ad in adapt:
        print(ad.__dict__)
