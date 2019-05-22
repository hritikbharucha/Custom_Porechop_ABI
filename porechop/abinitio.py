"""
Addition created by Quentin Bonenfant (quentin.bonenfant@gmail.com
                                       quentin.bonenfant@univ-lille.fr)
https://github.com/TheGTeck

This script is a dirty way to use my algorithm to find adapter sequence from the reads instead of
using the adapter.py static database. This is merely a test, but do seems to work fairly well.

A futur, cleaner version of this script will be edited in order to use C++ files "the proper way",
and maybe allowing users to specify arguments like k-mer size, k-mer threshold, sample size...
"""

import os
from  .adapters import Adapter

def  exportReadCheck(filename, check_reads):
    """
    EXPORT SEQUENCES IN CHECK_READ LIST TO FASTA
    (yeah... i know i know, don't hurt me please)
    """
    export_len = 100 #default length of ends to export
    starts = open(filename + "starts.fasta", "w")
    ends   =  open(filename +"ends.fasta" , "w")
    for read in check_reads:
        starts.write(">" + read.name + "\n")
        ends.write(">" + read.name + "\n")
        starts.write(read.seq[:export_len] + "\n")
        ends.write(read.seq[-export_len:] + "\n")
    starts.close()
    ends.close()

def execFindAdapt(check_reads):
    """ Manage the execution of compiled C++ files which
    are in charge of finding the adapter
    """

    adapter = {}

    # Temporary files prefix
    filename_pref = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tmp/temp_file_")
    tmpDir =  os.path.dirname(filename_pref) 
    # Creating tmp folder if not existing
    if( not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)
    # No need to clear, everything will be overwritten, or may lead to
    # funny results in the case of malfunction...align-content...  meh


    # exporting read check to fasta
    print("Processing read sample")
    exportReadCheck(filename_pref,check_reads)

    # Working on both ends
    for which_end in ['starts','ends']:

        # Counting kmers
        print("Counting k-mers from read sample " + which_end)
        os.system("./bin/kExtract " + filename_pref + which_end + ".fasta" + " -k 16 -lc 24 | head -n 500 > "+ filename_pref + which_end + "_kmers_first_500.fasta ")

        #Counting k-mers at ~ 2 errs
        print("Re-counting most frequent k-mers, allowing at most 2 err")
        os.system("./bin/adaptFinder3  " + filename_pref + which_end + ".fasta" + " -kf " + filename_pref + which_end + "_kmers_first_500.fasta -nt 4 -o " + filename_pref + which_end + "_COMPTAGE.txt")

        # Overlap step for reconstuction
        print('Reconstucting adapter from k-mers count')
        os.system("./bin/kOverlap " + filename_pref + which_end + "_COMPTAGE.txt")

        # extracting found adapter
        print("Fetching adapter sequence")
        with open( filename_pref + which_end + "_COMPTAGE_KOVERLAP.csv",'r') as f:
            # Trimming last 5 bases of adapter, which are often variable 
            adapter[which_end] = f.readline().split("\t")[0][:-5]

    # Getting back to porechop objects
    print("Building adapter object")
    adp = [Adapter("abinitio_adapter",
        start_sequence=('abinitio_Top' , adapter['starts']),
        end_sequence=('abinitio_Bottom', adapter['ends']))]

    print("Ab-initio phase: done")
    
    return(adp)