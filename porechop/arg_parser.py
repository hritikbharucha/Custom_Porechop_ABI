"""
Last modified 2022-03-22
Author: Quentin Bonenfant (quentin.bonenfant@gmail.com)
This file is a modified argument parser from Porechop.


Copyright 2017-2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop
This module contains the argument parser for Porechop.
This file is part of Porechop. Porechop is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop. If not, see <http://www.gnu.org/licenses/>.
Original Author: Ryan Wick
"""

import argparse
import multiprocessing
from .misc import MyHelpFormatter
from .version import __version__
import sys


def get_arguments():
    """
    Parse the command line arguments.
    """
    default_threads = min(multiprocessing.cpu_count(), 16)

    parser = argparse.ArgumentParser(description='Porechop_ABI: ab initio version of Porechop. '
                                                 'A tool for finding adapters in Oxford '
                                                 'Nanopore reads, trimming them from the ends and '
                                                 'splitting reads with internal adapters',
                                     formatter_class=MyHelpFormatter,
                                     add_help=False)

    abi_group = parser.add_argument_group('Ab-Initio options')
    abi_group.add_argument('-abi', '--ab_initio', action='store_true',
                           help='Try to infer the adapters from the read set '
                           'instead of just using the static database.')
    abi_group.add_argument('-go', '--guess_adapter_only', action='store_true',
                           help='Just display the inferred adapters then quit.')
    abi_group.add_argument('-abc', '--ab_initio_config', type=str,
                           help='Path to a custom config file for the '
                           'ab_initio phase (default file in Porechop folder)')
    abi_group.add_argument('-ws', '--window_size', type=int, default=3,
                           help='Size of the smoothing window used in the '
                           'drop cut algorithm. (default is 3, set to 1 to disable).')
    abi_group.add_argument('-ndc', '--no_drop_cut', action='store_true',
                           help='Disable the drop cut step entirely')
    abi_group.add_argument('-cap', '--custom_adapters', type=str,
                           help='Path to a custom adapter text file, '
                                'if you want to manually submit some.')
    abi_group.add_argument('-ddb', '--discard_database', action='store_true',
                           help='Ignore adapters from the Porechop database. '
                                'This option require either ab-initio (-abi) '
                                'or a custom adapter (-cap) to be set.')

    consensus_group = parser.add_argument_group('Consensus mode options')
    consensus_group.add_argument('-mr', '--multi_run', type=int, default=1,
                           help='Number of time the approximate count must be '
                                'performed to generate a consensus. '
                                'Each count file is exported separately.')
    consensus_group.add_argument('-cr', '--consensus_run', type=int, default=20,
                           help='If using multi-run option, set the number of '
                                'additional run performed if no stable consensus '
                                'is immediatly found.')
    consensus_group.add_argument('-ec', '--export_consensus', type=str,
                           help='Path to export the intermediate adapters found in '
                                'consensus mode.')
    consensus_group.add_argument('-aax', '--all_above_x', type=int, default=0,
                           help='Only select consensus sequences if they are made '
                                'using at least x percent of the total adapters. '
                                'Default is 0 (all consensus are kept).')
    consensus_group.add_argument('-box', '--best_of_x', type=int, default=0,
                           help='Only select the best x consensus sequences '
                                'from all consensus found.')


    graph_group = parser.add_argument_group('Graphs options')
    graph_group.add_argument('--export_graph', type=str,
                           help='Path to export the graphs used for assembly '
                                '(.graphml format), if you want to keep them')


    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                                 'trimmed reads will be printed to stdout)')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz'],
                            default='auto',
                            help='Output format for the reads - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = lots, '
                                 '3 = full - output will go to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')
    main_group.add_argument('-t', '--threads', type=int, default=default_threads,
                            help='Number of threads to use for adapter alignment')

    barcode_group = parser.add_argument_group('Barcode binning settings',
                                              'Control the binning of reads based on barcodes '
                                              '(i.e. barcode demultiplexing)')
    barcode_group.add_argument('-b', '--barcode_dir',
                               help='Reads will be binned based on their barcode and saved to '
                                    'separate files in this directory (incompatible with '
                                    '--output)')
    barcode_group.add_argument('--barcode_threshold', type=float, default=75.0,
                               help='A read must have at least this percent identity to a barcode '
                                    'to be binned')
    barcode_group.add_argument('--barcode_diff', type=float, default=5.0,
                               help="If the difference between a read's best barcode identity and "
                                    "its second-best barcode identity is less than this value, it "
                                    "will not be put in a barcode bin (to exclude cases which are "
                                    "too close to call)")
    barcode_group.add_argument('--require_two_barcodes', action='store_true',
                               help='Reads will only be put in barcode bins if they have a strong '
                                    'match for the barcode on both their start and end (default: '
                                    'a read can be binned with a match at its start or end)')
    barcode_group.add_argument('--untrimmed', action='store_true',
                               help='Bin reads but do not trim them (default: trim the reads)')
    barcode_group.add_argument('--discard_unassigned', action='store_true',
                               help='Discard unassigned reads (instead of creating a "none" bin)')

    adapter_search_group = parser.add_argument_group('Adapter search settings',
                                                     'Control how the program determines which '
                                                     'adapter sets are present')
    adapter_search_group.add_argument('--adapter_threshold', type=float, default=90.0,
                                      help='An adapter set has to have at least this percent '
                                           'identity to be labelled as present and trimmed off '
                                           '(0 to 100)')
    adapter_search_group.add_argument('--check_reads', type=int, default=10000,
                                      help='This many reads will be aligned to all possible '
                                           'adapters to determine which adapter sets are present')
    adapter_search_group.add_argument('--scoring_scheme', type=str, default='3,-6,-5,-2',
                                      help='Comma-delimited string of alignment scores: match, '
                                           'mismatch, gap open, gap extend')

    end_trim_group = parser.add_argument_group('End adapter settings',
                                               'Control the trimming of adapters from read ends')
    end_trim_group.add_argument('--end_size', type=int, default=150,
                                help='The number of base pairs at each end of the read which will '
                                     'be searched for adapter sequences')
    end_trim_group.add_argument('--min_trim_size', type=int, default=4,
                                help='Adapter alignments smaller than this will be ignored')
    end_trim_group.add_argument('--extra_end_trim', type=int, default=2,
                                help='This many additional bases will be removed next to adapters '
                                     'found at the ends of reads')
    end_trim_group.add_argument('--end_threshold', type=float, default=75.0,
                                help='Adapters at the ends of reads must have at least this '
                                     'percent identity to be removed (0 to 100)')

    middle_trim_group = parser.add_argument_group('Middle adapter settings',
                                                  'Control the splitting of read from middle '
                                                  'adapters')
    middle_trim_group.add_argument('--no_split', action='store_true',
                                   help='Skip splitting reads based on middle adapters '
                                        '(default: split reads when an adapter is found in the '
                                        'middle)')
    middle_trim_group.add_argument('--discard_middle', action='store_true',
                                   help='Reads with middle adapters will be discarded (default: '
                                        'reads with middle adapters are split) (required for '
                                        'reads to be used with Nanopolish, this option is on by '
                                        'default when outputting reads into barcode bins)')
    middle_trim_group.add_argument('--middle_threshold', type=float, default=90.0,
                                   help='Adapters in the middle of reads must have at least this '
                                        'percent identity to be found (0 to 100)')
    middle_trim_group.add_argument('--extra_middle_trim_good_side', type=int, default=10,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "good" side')
    middle_trim_group.add_argument('--extra_middle_trim_bad_side', type=int, default=100,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "bad" side')
    middle_trim_group.add_argument('--min_split_read_size', type=int, default=1000,
                                   help='Post-split read pieces smaller than this many base pairs '
                                        'will not be outputted')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

    args = parser.parse_args()

    # Checking Scoring Scheme
    try:
        scoring_scheme = [int(x) for x in args.scoring_scheme.split(',')]
    except ValueError:
        sys.exit('Error: incorrectly formatted scoring scheme')
    if len(scoring_scheme) != 4:
        sys.exit('Error: incorrectly formatted scoring scheme, expected 4 values.')
    args.scoring_scheme_vals = scoring_scheme

    if args.barcode_dir is not None and args.output is not None:
        sys.exit(
            'Error: only one of the following options may be used: --output, --barcode_dir')

    if args.untrimmed and args.barcode_dir is None:
        sys.exit('Error: --untrimmed can only be used with --barcode_dir')

    if args.barcode_dir is not None:
        args.discard_middle = True

    if args.output is None and args.barcode_dir is None:
        args.print_dest = sys.stderr
    else:
        args.print_dest = sys.stdout

    if args.threads < 1:
        sys.exit('Error: at least one thread required')

    # Checking if number of ab_initio run is above 1
    if(args.multi_run < 1):
        print("Ab-Initio multi-runs can only be set at 1 or above.")
        print("Setting to 1")
        args.multi_run = 1
    # Same for consensus runs
    if(args.consensus_run < 1):
        print("Ab-Initio supplementary (consensus) runs can not be below 1")
        print("Setting to 1")
        args.consensus_run = 1

    # Force setting ab_initio if we only want to find adapter.
    if(args.guess_adapter_only):
        args.ab_initio = True

    # Checking if at least one adapter is in the database
    if(args.discard_database):
        if(not(args.ab_initio or args.custom_adapters)):
            print("At least one adaptater is needed to make Porechop run")
            print("Please either set -abi option for ab-initio inferrence or")
            print("specifie a custom adapter set using -cap and a file.")
            exit(1)
    return args
