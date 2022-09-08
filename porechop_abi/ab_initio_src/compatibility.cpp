#include <iostream>
#include "compatibility.h"



//////////////////////////////////////////////////////////////////////////////
// Semi-global alignment of pre-loaded align object with free gaps at the end
void semi_global_align(TAlign & align){
    // score is not used, so we don't return it
    globalAlignment(align, Score<int, Simple>(match, mismatch, gap), AlignConfig<true, true, true, true>());
}


//////////////////////////////////////////////////////////////////////////////
// Build SeqAn DNA String Set from two sequences, longest sequence first
TStringSet make_stringSet(TSequence seq1, TSequence seq2){
    // Storing in string set, longest sequences first.
    TStringSet sequences;
    if(length(seq1) < length(seq2)){
        appendValue(sequences, seq2);
        appendValue(sequences, seq1);
    }
    else{
        appendValue(sequences, seq1);
        appendValue(sequences, seq2);
    }
    return(sequences);
}



//////////////////////////////////////////////////////////////////////////////
// Finding mapping start and end.
// This function is inspired by Porechop alignment.cpp.
void get_alignment_parameters(TAlign & align, int * parameters){
    // Row access:
    TRow row1 =  row(align, 0);
    TRow row2 =  row(align, 1);

    int alignmentLength = length(row1);

    int alignmentStartPos = -1;
    bool seq1_started = false;
    bool seq2_started = false;
    for (int i = 0; i < alignmentLength; ++i) {
        if (getValue(row1, i) != '-')
            seq1_started = true;
        if (getValue(row2, i) != '-')
            seq2_started = true;
        if (seq1_started && seq2_started) {
            alignmentStartPos = i;
            break;
        }
    }
    parameters[0] = alignmentStartPos;

    // We use the same logic to see when the alignment has ended.
    int alignmentEndPos = -1;
    bool readEnded = false;
    bool adapterEnded = false;
    for (int i = alignmentLength - 1; i >= 0; --i) {
        if (getValue(row1, i) != '-')
            readEnded = true;
        if (getValue(row2, i) != '-')
            adapterEnded = true;
        if (readEnded && adapterEnded) {
            alignmentEndPos = i;
            break;
        }
    }
     parameters[1] = alignmentEndPos;
}


//////////////////////////////////////////////////////////////////////////////
// Get the edit distance between two sequences on the ALIGNED area only.
int get_edit_distance(TAlign & align, int alignmentStartPos, int alignmentEndPos){
    // getting row handles back
    TRow row1 = row(align, 0);
    TRow row2 = row(align, 1);

    // Computing distance
    int levenstein_distance = 0;
    for(int i = alignmentStartPos; i < alignmentEndPos; i++){
        if(getValue(row1, i) != getValue(row2, i))  
            levenstein_distance ++;
    }
    return(levenstein_distance);
}


//////////////////////////////////////////////////////////////////////////////
// Determine if seq2 is included in seq1.
bool is_included(TAlign & align, int alignmentStartPos, int alignmentEndPos){

    // getting row handles back
    TRow row1 = row(align, 0);
    TRow row2 = row(align, 1);
    // getting StringSet handle
    TStringSet sequences = stringSet(align);

    // seq2 is always smaller than seq1 so we only have to check if seq2 is in seq1
    int s1_start = toSourcePosition(row1, alignmentStartPos);
    int s1_end = toSourcePosition(row1, alignmentEndPos);
    int s1_len = length(sequences[0]);
    // int s2_start = toSourcePosition(row2, alignmentStartPos);
    // int s2_end = toSourcePosition(row2, alignmentEndPos);
    // int s2_len = length(sequences[1]);

    // if either start position on the alignment or end position on the alignment are not close to the ends of s1,
    // we consider the sequence 2 as "inside" sequence 1.
    return(( s1_start > OVERLAP_LIMIT  or (s1_len - s1_end) > OVERLAP_LIMIT));
}


//////////////////////////////////////////////////////////////////////////////
// Return the compatibility flag of two aligned sequences
// Output flags depends on the alignement state:
// 0: not compatible (score too high)
// 1: compatible
// 2: compatible, but s2 is included in s1
int get_flag(TAlign & align){

    // Alignment  (score discarded here)
    semi_global_align(align);

    // Finding alignment parameters
    int align_param [2];
    get_alignment_parameters(align, align_param);
    int alignmentStartPos = align_param[0];
    int alignmentEndPos = align_param[1];

    // Get mapped length from alignement start and end
    int mappedLength = alignmentEndPos - alignmentStartPos + 1;

    // Counting errors between the alignment boundaries
    int levenstein_distance = get_edit_distance(align, alignmentStartPos, alignmentEndPos);

    // Detecting if one sequence is included in the other:
    bool included = is_included(align, alignmentStartPos, alignmentEndPos);

    // Identity (in percent) of the aligned region.
    float identity = (mappedLength - levenstein_distance) * 100 / mappedLength;

    // Finally, define the Flag, default is: not compatible
    int flag = 0;
    // Identity flag
    if(identity >= identity_threshold){
        flag = 1;
    }
    // Inclusion flag
    if(included and flag == 1){
        flag = 2;
    }
    return(flag);
}


// Entry point function to get the compatibility flag when using this as an external library.
int check_compatibility(TRawSeq raw_seq1, TRawSeq raw_seq2){
    // Sequence building
    TSequence seq1 = raw_seq1;
    TSequence seq2 = raw_seq2;
    TStringSet sequences;
    sequences = make_stringSet(seq1, seq2);
    TAlign align(sequences);
    int flag = get_flag(align);
    return(flag);
}