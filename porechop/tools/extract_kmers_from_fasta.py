"""
Pase a fasta file containing sequences and extract the k-mers from
the start / end, then output them to stdout.
Author: Quentin Bonenfant (quentin.bonenfant@gmail.com)
"""
import sys
import os
import re

DNA_rule = re.compile(r"^[ATCG]+$")
RNA_rule = re.compile(r"^[AUCG]+$")

def print_usage():
    print("Usage: extract_kmers_from_fasta <fasta> <k> [s] [e] > outfile")
    print("fasta: the path to the fasta file")
    print("k: the k-mers size")
    print("Optional:")
    print("s : number of bases from the start of the sequence")
    print("e : number of bases from the end of the sequence")
    print("Note:")
    print("s and e must be given together or not at all")
    print("If not set, all the k-mers are printed")
    exit(1)


def get_kmers(seq, k):
    """Return the set of k-mers from a sequence
    @param seq: the sequence (string)
    @param k : the k-mer size
    @return kms: The k-mers set (set)
    """
    kms = set()
    for i in range(len(seq) - k + 1):
        kms.add(seq[i:i+k])
    return(kms)


def is_nucleotides(seq):
    """Predicate checking if seq is a DNA or RNA sequence
    @param seq: a potential DNA / RNA sequence (string)
    @return Is the sequence DNA / RNA ? (bool)
    """
    return(DNA_rule.match(seq) is not None or RNA_rule.match(seq) is not None)


def fasta_reader(filename):
    # checking if filename is a file
    if(not os.path.isfile(filename)):
        print("Input is not a Fasta file.")
        print_usage()
    # If ok, parse the file
    with open(filename) as f:
        seq = ""
        for line in f:
            if(line[0] == ">"):
                if(seq):
                    if(is_nucleotides(seq)):
                        yield(seq)
                    else:
                        print("/!\\\tA sequence in the fasta file is not",
                              file=sys.stderr)
                        print("a DNA or RNA sequence.",
                              file=sys.stderr)

                        print("-> " + seq[:40] + "...",
                              file=sys.stderr)
                        exit(1)
                seq = ""
            else:
                seq += line.rstrip("\n")
        if(is_nucleotides(seq)):
            yield(seq)
        else:
            print("/!\\\tA sequence in the fasta file is not",
                  file=sys.stderr)
            print("a DNA or RNA sequence.",
                  file=sys.stderr)
            print("-> " + seq[:40] + "...",
                  file=sys.stderr)
            exit(1)


if __name__ == '__main__':
    nb_args = len(sys.argv)
    if( nb_args < 3):
        print_usage()
    try:
        k = int(sys.argv[2])
        infile = sys.argv[1]
    except ValueError:
        print("k-mer size need to be an int.")
        print_usage()
    start = 0
    end = 0
    need_split = False

    # If start / end limit are specified
    if(nb_args == 5):
        try:
            start =int(sys.argv[3])
            end = int(sys.argv[4])
        except ValueError:
            print("Could not process start / end limits (int needed).",
                  file=sys.stderr)
            print("Using the whole sequence instead.",
                  file=sys.stderr)
        else:
            need_split = True

    kms = set()
    for seq in fasta_reader(infile):
        if(need_split):
            kms |= get_kmers(seq[:start], k)
            kms |= get_kmers(seq[-end:], k)
        else:
            kms |= get_kmers(seq, k)
    for km in kms:
        print(km)
