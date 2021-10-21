# Approx_counter

This is a tool designed to count k-mers from the ends of the sequences of a multi fasta/fastq file, allowing for up to 2 errors.
It was developped to be paired with Porechop_ABI, to help infer the sequences of [Oxford Nanopore](https://nanoporetech.com/) adaptaters using approximate kmer count.
Since adapter are short sequence at the ends of the read, the kmer composing them should be extremely frequent in those area.
That's why we count kmers allowing an edit distance up to 2 between the ressearched kmer and an indexed sample of the reads.
My goal here is to provide a simple tool able to perform this task.

A simple assembly method can then be used to reconstruct the potential adapter, more information can be found in the Porechop_ABI repo.


## Requirement
- C++ (11+)
- SeqAn (2.4.0+)

## Compiling
In order to compile this file, I recommend using the following command
~~~
g++ -std=c++14 -fopenmp  -O3 -DNDEBUG -march=native  -mtune=native  approx_counter.cpp -lrt -o approx_counter
~~~

## Usage
REQUIRED ARGUMENTS

    input_filename STRING

OPTIONS

    -h, --help
          Display the help message.
    -lc, --low_complexity DOUBLE
          low complexity filter threshold (for k=16), default 1.5
    -sn, --sample_n INTEGER
          sample n sequences from dataset, default 10k sequences
    -sl, --sample_length INTEGER
          size of the sampled portion, default 100 bases
    -nt, --nb_thread INTEGER
          Number of thread to work with, default is 4
    -k, --kmer_size INTEGER
          Size of the kmers, default is 16
    -lim, --limit INTEGER
          limit the number of kmer used after initial counting, default is 500
    -v, --verbosity INTEGER
          Level of details printed out (fixed for the moment)
    -e, --exact_file STRING
          path to export the exact k-mer count, if needed. Default: no export
    -o, --out_file STRING
          path to the output file, default is ./out.txt

## Example
approx_counter file.fasta -k 16 --sample_n 20000 --sample_length 90 -nt 4 lim 1000 -e exact_out.txt -o approx_out.txt


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)


