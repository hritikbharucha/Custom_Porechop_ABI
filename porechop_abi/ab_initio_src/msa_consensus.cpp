#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>  

using namespace seqan;

// Default scoring parameters (same as Porechop)
#define match 5
#define mismatch -3
#define gap -1
#define gapExt -3


int main(int argc, char const ** argv)
{   /////////////////////////////////////////////////
    // Parsing the input fasta file
    //
    if(argc < 2){
        std::cerr << "Usage: msa_consensus <fasta file>" << std::endl;
        return(1);
    }
    std::string input_file = argv[1];
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;
    SeqFileIn seqFileIn(toCString(input_file));
    readRecords(ids, seqs, seqFileIn);
    int nb_seq = length(seqs);


    /////////////////////////////////////////////////
    // Building alignment structures
    //
    Align<DnaString> align;
    resize(rows(align), nb_seq);
    for (int i = 0; i < nb_seq; ++i)
        assignSource(row(align, i), seqs[i]);

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

    /////////////////////////////////////////////////
    // Create the profile string
    String<ProfileChar<Dna> > profile;
    resize(profile, length(row(align, 0)));
    for (int rowNo = 0; rowNo < nb_seq; ++rowNo)
        for (unsigned i = 0; i < length(row(align, rowNo)); ++i)
            profile[i].count[ordValue(getValue(row(align, rowNo), i))] += 1;

    ////////////////////////////////////////////////
    // Call consensus from this string
    DnaString consensus;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        int idx = getMaxIndex(profile[i]);
        if (idx < 4)  // is not gap
            appendValue(consensus, Dna(getMaxIndex(profile[i])));
    }
    ///////////////////////////////////////////////
    // Export consensus sequence
    std::cout << consensus <<std::endl;

    return 0;
}