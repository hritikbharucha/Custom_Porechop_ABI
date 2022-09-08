#ifndef COMPATIBILITY_H
#define COMPATIBILITY_H

#include <string>
#include <seqan/align.h>
#include <seqan/sequence.h>

// Default scoring parameters
#define match 2
#define mismatch -1
#define gap -1
#define gapExt -1

// minimum identity to build a link.
#define identity_threshold 87.5

// Number of bases required to considere a sequence as 'part' of the other
#define OVERLAP_LIMIT 3

using namespace seqan;

// Types definition
typedef int Tscore;
typedef String<Dna> TSequence;               // sequence type
typedef StringSet<CharString> TCharSet;      // Type for CharStrng string sets. 
typedef StringSet<TSequence> TStringSet;     // container for strings
typedef Align<TSequence, ArrayGaps> TAlign;  // align type for semi-global (AlignGap based alignement method)
typedef Row<TAlign>::Type TRow;              // gapped sequence type
typedef char *  TRawSeq;                     // Raw sequences from Python
typedef const char *  TRawConstSeq;          // Raw sequences from stdin



// Semi-global alignment of pre-loaded align object with free gaps at the end
void semi_global_align(TAlign & align);

// Convert the raw input from stdin to SeqAn DNA String Set.
TStringSet make_stringSet(TSequence seq1, TSequence seq2);

// Build alignments directly from characters
TAlign build_align(TRawSeq raw_seq1, TRawSeq raw_seq2);

// Finding mapping start and end.
void get_alignment_parameters(TAlign & align, int * parameters);

// Get the edit distance between two sequences on the ALIGNED area only.
int get_edit_distance(TAlign & align, int alignmentStartPos, int alignmentEndPos);
    
// Determine if seq2 is included in seq1.
bool is_included(TAlign & align, int alignmentStartPos, int alignmentEndPos);

// Return the compatibility flag of two aligned sequences
int get_flag(TAlign & align);

// External access for Python script, need C linkage.
extern "C" {
    int check_compatibility(TRawSeq raw_seq1, TRawSeq raw_seq2);
}
#endif // COMPATIBILITY_H