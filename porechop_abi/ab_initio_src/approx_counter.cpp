#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>  
#include <seqan/arg_parse.h>

#include <iostream>
#include <cstdlib>  
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <stdexcept>
#include <unordered_map>
#include <set>

using namespace seqan;

// Program start timestamp
const auto boot_time = std::chrono::steady_clock::now();

// Alphabet used for 2 bit conversion
const std::string DNA = "ACGT";

// Max number of errors, need to be fixed at compile time
const uint8_t MAXERR = 2; 


// Setting the index
typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
using index_t = Index<StringSet<Dna5String>, BidirectionalIndex<FMIndex<void,TFastConfig> > >;

// counter type, using unordered map.
using counter = std::unordered_map<uint64_t,uint64_t>;
// pair vector
using int_pair = std::pair<uint64_t,uint64_t>;
using pair_vector = std::vector<int_pair> ;
// type of sequence set
using sequence_set_type = StringSet<Dna5String> ;
// vector of boolean used to keep track of kmer positions count.
using bit_field = std::vector<bool>  ;
// config file parameter map
using arg_map = std::unordered_map<std::string, std::string>;
// Set of kmers (2bit representation)
using kmer_set_t = std::set<uint64_t>;


/**
    Convert a Seqan Dna5String to uint64_t int.
    SeqAn already store Dna5String in 2 bit representation,
    but it is not easy to access it "as is" or to use it
    as key in a hash map.
    @param The sequence to convert
    @return the kmer in 2 bit format, as an uint64_t int.
*/
inline uint64_t dna2int(Dna5String seq){
    
    uint64_t value = 0;
    for(auto c : seq){
        value = value << 2 | (uint8_t)(c);
    }
    return(value);
}

/**
    Convert an uint64_t int back to Seqan Dna5String
    @param the integer to convert
    @param the size of the kmer
    @return the kmer in SeqAn Dna5String format
*/
inline Dna5String int2dna(uint64_t value, uint8_t k){
    std::string seq = "";
    for(int i = 0; i< k; i++){
        seq = DNA[value & 3] + seq;
        value >>= 2;

    }
    return(Dna5String(seq));
}


/**
    Shortcut to print text in stdout with a time stamp.
    @param whatever you want to print.
*/
template<typename TPrintType>
void print(TPrintType text, int tab = 0)
{
    const auto milis = std::chrono::duration <double, std::milli>(std::chrono::steady_clock::now() - boot_time).count();
    std::cout << "[" << milis << " ms]\t" ;
    for(int i = 0; i < tab; i++){
        std::cout << "\t";
    }
    std::cout << text << std::endl;
}

/**
    Extremly simple config file parser.
    format  : args=value , one per line
    comments: #  on line start
    @param the path to config file
    @return a map containing the set parameters.
*/
arg_map parse_config(std::string inputFile){

    arg_map params;

    std::ifstream inFile;
    inFile.open(inputFile);
    if(inFile.is_open()){
        std::string line;
        while( std::getline(inFile, line) ){
            std::string arg = "";
            std::string val = "";
            bool sep = false;
            if(line[0] != '#'){
                for(auto c: line){
                    if(c == '='){
                        sep = true;
                    }
                    else if( c!= ' '){
                        if(sep)
                            val+=c;
                        else
                            arg+=c;
                    }
                }
                params[arg] = val;
            }
        }
    }
    else{
        std::cerr << "/!\\ WARNING: Could not open config file\n";
    }
    return(params);
}


/**
    Shortcut to print list of kmer counts (using either pairs or map )
    @param the counter to export (pair or map of kmer associated to their count)
    @param k, the size of the kmer counted (needed for conversion 2bit representation back to DNA)
*/
template<typename TPrintType>
void printCounters(TPrintType & pvec, uint8_t k){
    for(auto it = pvec.begin(); it!= pvec.end(); ++it)
    { 
        std::cout << int2dna(it->first,k) << " " << it->second << "\n";
    }   
}

/**
    Export a counter to a file
    @param the counter to export (pair or map of kmer associated to their count)
    @param k, the size of the kmer counted (needed for conversion 2bit representation back to DNA)
    @param the path to the outputfile
*/
template<typename TPrintType>
bool exportCounter(TPrintType & pvec, uint8_t k, std::string output){

    std::ofstream outputFile;
    outputFile.open (output);
    if(outputFile.is_open()){
        for(auto it = pvec.begin(); it!= pvec.end(); ++it)
        { 
            outputFile << int2dna(it->first,k) << "\t" << it->second << "\n";
        }   
        outputFile.close();
    }
    else{
        std::cerr << "/!\\ ERROR: COULD NOT OPEN FILE " << output << std::endl;
        return(false);
    }
    return(true);
}

/**
    Adjust low complexity threshold value to kmer size
    @param low complexity threshold for a size of kmer
    @param kmer size for this threshold
    @param new k-mer size
    @return the appropriate threshold value for the new kmer size.
*/
float adjust_threshold(float c_old, uint8_t k_old, uint8_t k_new ){
    float c_new = c_old * float( std::pow(k_new - 2 + 1,2) /  std::pow(k_old - 2 + 1,2));
    return(c_new);
}

/**
    Compute the sum of a vector (of int for example).
    @param the vector
    @return the sum
*/
template<typename TVector>
inline uint64_t vectorSum(TVector vec)
{
    uint64_t res = 0;
    for(auto it: vec)
    {
        res += it;
    }
    return(res);
}

/**
    Check the complexity of a kmer
    Our score is derived from 2006 DUST publication on
    fast DNA sequence masking
    https://doi.org/10.1089/cmb.2006.13.1028
    @param the kmer to test, in 2 bit representation, cast as an uint64_t int
    @param k the size of said kmer
    @param threshold for the low complexity filter.
    @return True if the kmer contains low complexity region
*/
inline  bool haveLowComplexity(uint64_t kmer, uint8_t k, float threshold){
    
    uint64_t counts[16] = { 0 }; // 16 possibles dimers
    // reading using sliding window of 2
    for(int i = 0; i < k-1; i++){
        // storing dimers as 2 * 2 bits
        uint8_t c =  kmer & 15;
        // removing last element of the k-mer
        kmer >>=2;
        // updating value of dimer count
        counts[c]++;
    }

    float s = 0;
    size_t sum = 0;
    for(auto v:counts){
        sum +=  v * (v-1);  
    }
    s =  sum / float(2 * (k-2));
    return s>= threshold;
}



/**
    Return the complexity of a kmer
    Our score is derived from 2006 DUST publication on
    fast DNA sequence masking
    https://doi.org/10.1089/cmb.2006.13.1028
    @param the kmer to test, in 2 bit representation, cast as an uint64_t int
    @param k the size of said kmer
    @return True if the kmer contains low complexity region
*/
float getComplexity(uint64_t kmer, uint8_t k){
    
    uint64_t counts[16] = { 0 }; // 16 possibles dimers
    // reading using sliding window of 2
    for(int i = 0; i < k-1; i++){
        // storing dimers as 2 * 2 bits
        uint8_t c =  kmer & 15;
        // removing last element of the k-mer
        kmer >>=2;
        // updating value of dimer count
        counts[c]++;
    }

    float s = 0;
    size_t sum = 0;
    for(auto v:counts){
        sum +=  v * (v-1);  
    }
    s =  sum / float(2 * (k-2));
    return s;
}

/**
    Object used to compare two k-mer count.
    Such object is needed because k-mer size is
    a required parameter that can't just be passed
    to a delegate function.
*/
struct CompareCount{
    // Constructor
    CompareCount(int k) { this->k = k; }    
    /*  Kmer counter comparison function
    @param a, first count object
    @param b, second count object
    @return is a smaller than b
    */
    bool operator () (int_pair a, int_pair b){
        // If counts are equals
        if(a.second == b.second){
            // sort by complexity
            float a_comp = getComplexity(a.first, this->k);
            float b_comp = getComplexity(b.first, this->k);
            // if still equal, use decreasing lexicographic order
            if(a_comp == b_comp){
                return(a.first > b.first);
            }
            else{
                return(a_comp < b_comp);    
            }
            
        }
        // if not, just regular count-based comparison.
        else{
            return(a.second>b.second);
        }
    }
    // Defining local k
    int k;
};


/**
    Check if a DNA5 sequence only contains ATCG
    @param the DNA5 sequence
    @return True if the sequence do not contains IUPAC chars (only "ACGT" chars).
 */
bool is_DNA(Dna5String seq){
    for(auto c: seq){
        // If the c is a "N", int representation will be >= 4.
        if(uint8_t(c) >= 4){
            return(false);
        }
    }
    return(true);
}


/**
    Check if the kmer is autorised or not
    @param the kmer to test, in 2 bit representation, cast as an uint64_t int
    @param a set containing the forbidden kmers
    @return True if the kmer is contained in the set
*/
inline  bool isForbiddenKmer(uint64_t kmer, kmer_set_t & kmer_set){
    return( kmer_set.find(kmer) != kmer_set.end() ) ;
}

/**
    Parse a kmer list file and return a kmer set.
    Only true DNA kmers are accepted (no IUPAC).
    @param the file containging a list of kmer to exclude
    @return A kmer set (set of uint64_t)
*/
kmer_set_t parse_kmer_list(std::string kmer_file){
    
    kmer_set_t kmer_set;

    std::ifstream kmer_file_stream(kmer_file);
    // parsing file
    if( kmer_file_stream.is_open()) {
        for( std::string line; getline( kmer_file_stream, line );){
            
            // Converting line to DNA5, invalid chars will be converted to N.
            Dna5String km = line;

            // Inserting new kmer only if it is standard DNA (ATCG)
            if(is_DNA(km)){
                kmer_set.insert(dna2int(km));
            }
        }
        kmer_file_stream.close();
    }
    else{
        std::cerr << "/!\\ ERROR: COULD NOT OPEN EXCLUDED KMER FILE, must quit\n";
        exit(1);
    }
    return(kmer_set);
}

/**
    Return the solid kmers. Kmer with count over threshold will be kept.
    @param kmer count map (see count_kmer and error count)
    @param The minimum abundance of the kmers to return
    @return vector of pair containing the solids kmers and their associated count.
*/
pair_vector get_solid_kmers(counter & count_map, uint64_t solid_km){

    pair_vector kmer_vec(std::make_move_iterator(count_map.begin()), std::make_move_iterator(count_map.end()));
    std::sort(kmer_vec.begin(), kmer_vec.end(), [](auto x, auto y){ return x.second > y.second;} );
    uint64_t limit = 0;
    // Finding the position in the vector where kmer count is no longer high enough
    for(auto kmer_count: kmer_vec){
        if(kmer_count.second >= solid_km ){
            limit ++;
        }
        else{
            break;
        }
    }
    kmer_vec.resize(limit);
    return(kmer_vec);
}

/**
    Return the first top kmers, ranked by count, descending.
    @param kmer count map (see count_kmer and error count)
    @param number of kmers to return
    @return vector of pair containing the most frequent kmers and their associated count.
*/
pair_vector get_most_frequent(counter & count_map, uint64_t limit, int k){

    pair_vector kmer_vec(std::make_move_iterator(count_map.begin()), std::make_move_iterator(count_map.end()));
    //std::sort(kmer_vec.begin(), kmer_vec.end(), [](auto x, auto y){ return x.second > y.second;} );
    std::sort(kmer_vec.begin(), kmer_vec.end(), CompareCount(k) );
    if(kmer_vec.size() > limit){
        kmer_vec.resize(limit);
    }
    return(kmer_vec);
}


/**
    Sample the sequences set and return requested samples cut to size
    @param Set of sequences to sample (SeqAn StringSet of Dna5String)
    @param Number of sequences to sample
    @param Size of the sampled sequence.
    @return A set of sample sequences cut to size.
*/
sequence_set_type sampleSequences(sequence_set_type & sequence_set, uint64_t nb_sample, uint64_t cut_size, bool bot, uint8_t v){
    sequence_set_type sample;
    
    // Initialising the random seed
    srand(time(0));

    uint64_t sequence_set_size = length(sequence_set);
    // Building vector with all possible seq indice
    std::vector<int> vec(sequence_set_size) ; 
    std::iota(std::begin(vec), std::end(vec), 0);

    // and applying random shuffling to said vector.
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(vec.begin(), vec.end(), g);

    // counters
    uint64_t nb_seq = 0;
    uint64_t i = 0;
    uint64_t seq_id;

    // display
    if(v>0){
        if(bot){
            print("Sampling the ends of reads", 1);
        }
        else{
            print("Sampling the start of reads", 1);
        }
    }

    // Fetching the random sequences
    while(nb_seq < nb_sample and i < sequence_set_size ){
        
        // current sequence id
        seq_id = vec[i];
        
        // asjusting cut size to read length, if it is too short.
        uint64_t current_cut_size = std::min(uint64_t(length(sequence_set[ seq_id ])), cut_size);

        if(current_cut_size < cut_size and v >= 2){
            std::cerr << "/!\\ WARNING: Cut size is longer that current read! (read id: " << seq_id << ")."<< std::endl;
        }

        // we also need to check that the sequence is
        // at least long enough to contains both adapters.
        if( length(sequence_set[ seq_id ]) >= cut_size * 2 ){
            if(bot){
                appendValue(sample, suffix(sequence_set[ seq_id ], length(sequence_set[ seq_id ]) - 1 - current_cut_size ));
            }
            else{
                appendValue(sample, prefix(sequence_set[ seq_id ], current_cut_size));
            }
            
            nb_seq +=1;
        }
        i++;
    }
    if(v>0)
        print("Sampled " + std::to_string(length(sample)) + " sequences", 1);
    return(sample);
}


/**
    Perform a simple exact count of all the k-mers from a sample set of sequences
    @param Set of sequences (SeqAn StringSet of Dna5String)
    @param k, size of the kmers
    @param threshold of the low complexity filter
    @return a map of the kmer count, with a kmer hash as key.

*/
counter count_kmers(sequence_set_type & sequences, uint8_t k, float threshold, kmer_set_t & kmer_set){

    counter count;
    unsigned had_n = 0; // keeping track of k-mers with N
    for(auto seq: sequences){
        unsigned i = 0;
        uint64_t n; // storing k-mers in 2 bit format as 64 bit int using C++ std types.
        // First kmer
        Infix<Dna5String>::Type km = infix(seq, i, i+k);
        while(i+k <= length(seq)){ 
            // If current k-mer does not contains N or other IUPAC chars, proceed
            if(is_DNA(km)){
                n = dna2int(km); // converting k-mer to int representation
                // If the k-mers is not forbidden or low complexity, add to count.
                if(not haveLowComplexity(n, k,threshold) and not isForbiddenKmer(n, kmer_set)){
                    count[n] +=1;
                }
            }
            else{
                had_n ++;
            }
            // Updating k-mer
            i++;
            km = infix(seq, i, i+k);
        }
    }
    if(had_n > 0){
        std::cerr << "/!\\ WARNING: This dataset contained sequences with 'N' symbols. ";
        std::cerr << "/!\\ WARNING: Current implementation ignores k-mers containing 'N'.";
        std::cerr << "/!\\ WARNING: A total of "<< had_n << " k-mers were ignored." << std::endl;
    }
    return(count);
}


/**
    Search and count a list of kmer in a set of sequences, at a Levenstein distance of at most 2.
    @param Set of sequences (SeqAn StringSet of Dna5String)
    @param the previous count of exact kmer (the kmer list)
    @param number of thread to use
    @param k, size of the kmers
    @return a map of the kmer count, with a kmer hash as key.

*/
counter errorCount( sequence_set_type & sequences, pair_vector & exact_count, uint8_t nb_thread, uint8_t k, uint8_t v){
    

    uint64_t sample_size = length(sequences);
    if(v>0)
        print("Preparing index",1);
    index_t  index(sequences);
    
    if(v>0)
        print("Creating index",1);
    indexCreate(index);
    
    // Result storage
    counter results;

    // setting number of parallel thread
    omp_set_num_threads(nb_thread);
    if(v>0)
        print("Starting approximate counting",1);
    #pragma omp parallel shared(index, results)
    {
        // local variable to keep track of kmer occurences
        std::array<bit_field,3> tcount;

        // Delegate function for SeqAn find function (process occurences)
        auto delegateParallel = [& tcount](auto & iter, const Dna5String & needle, int errors)
        {
            (void)needle; // casting to void to avoid "unused variable" warning
            // Needle is required for the delegate to be accepted but is not needed in this program.

            for (auto occ : getOccurrences(iter)){
                uint64_t read_id = getValueI1(occ);
                tcount[errors][read_id] = true;
                }
        };

        #pragma omp for schedule(dynamic)
        for(unsigned km_id=0; km_id<length(exact_count); km_id++)
        {
            
            /*  Supposing a kmer is very unlikely to be twice in the same read start
                we just store the read id in which kmer have been found using a bit field
                which size is the number of reads.
                Since we search at 2 error, we can either:
                    - Find the kmer exactly
                    - Find it at 1 error
                    - Find it at 2 errors
                So there is one bit field per number of error.
            */
            for(int i=0; i<3; i++){
                tcount[i] = bit_field(sample_size,false);
            }
            // 2 bit encoded kmer as uint64_t int
            uint64_t kmer = exact_count[km_id].first;
            // ressearch, filling tcount
            find<0, MAXERR >(delegateParallel, index, int2dna(kmer,k), EditDistance() );

            
            // computing total number of occurences
            uint64_t total = 0;
            for(auto bit_count: tcount){
                total +=  vectorSum(bit_count);
            }
            // Updating global counter
            #pragma omp critical
            results[kmer] = total;
                
        }
    }
    return(results);
}


seqan::ArgumentParser::ParseResult get_args(seqan::ArgumentParser & parser, int argc, char const ** argv){

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "input filename"));

    addOption(parser, seqan::ArgParseOption(
        "lc", "low_complexity", "low complexity filter threshold (for k=16), default 1.5",
        seqan::ArgParseArgument::DOUBLE, "kmer filename"));
    
    addOption(parser, seqan::ArgParseOption(
        "sn", "sample_n", "sample n sequences from dataset, default 10k sequences",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "sl", "sample_length", "size of the sampled portion, default 100 bases",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "nt", "nb_thread", "Number of thread to work with, default is 4",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_size", "Size of the kmers, default is 16",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "lim", "limit", "limit the number of kmer used after initial counting, default is 500",
        seqan::ArgParseArgument::INTEGER, "INT"));
        addOption(parser, seqan::ArgParseOption(
        "mr", "multi_run", "Number of time the count must be performed. Each count is exported separately.",
        seqan::ArgParseArgument::INTEGER, "INT"));

     addOption(parser, seqan::ArgParseOption(
        "v", "verbosity", "Level of details printed out",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "e", "exact_file", "path to export the exact k-mer count, if needed. Default: no export",
        seqan::ArgParseArgument::STRING, "exact count output file"));

    addOption(parser, seqan::ArgParseOption(
        "conf", "config", "path to the config file",
        seqan::ArgParseArgument::STRING, "config file"));

    addOption(parser, seqan::ArgParseOption(
        "fk", "forbidden_kmer", "take a file containing 'forbidden' kmers, excluding them from the search pool. One kmer per line.",
        seqan::ArgParseArgument::STRING, "config file"));

    addOption(parser, seqan::ArgParseOption(
        "sk", "solid_km", "Use solid kmer instead of most frequents. This option will override sample number (-sn / --sample_n).",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "se", "skip_end", "Skip end adapter ressearch (only search start). /!\\ If this option is set, and adaptFinder is run trough PorechopABI, the --guess_only / -go MUST be set."
        ));

    addOption(parser, seqan::ArgParseOption(
        "o", "out_file", "path to the output file, default is ./out.txt",
        seqan::ArgParseArgument::STRING, "output file"));

    // Parser command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    return(res);

}


///////////////////////////////////////////////////////////////////////////////
//
//          MAIN FUNCTION
//
///////////////////////////////////////////////////////////////////////////////


int main(int argc, char const ** argv)
{
    ////////////////////////////////////////////
    //
    //    ARGUMENT PARSER
    //
    ////////////////////////////////////////////

    ////////////////////////
    // SeqAn ArgumentParser
    ////////////////////////

    seqan::ArgumentParser parser("adaptFinder");
    // Parsing
    seqan::ArgumentParser::ParseResult res = get_args(parser, argc, argv);

    // If parsing was not successful then exit with code 1. if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Default values for parameters
    std::string output = "out.txt";     // output file
    std::string exact_out;   // exact count output file
    std::string config_file; // configuration file
    std::string forbid_kmer; // forbidden kmers file, one kmer per line.
    uint64_t solid_km= 0;    // Use solid k-mer instead of most frequent
    uint64_t nb_thread = 4;  // default number of thread
    uint64_t k = 16;         // kmer size, 2<= k <= 32
    uint64_t sl = 100 ;      // sequence sampling size
    uint64_t sn = 40000;     // number of sequence sampled
    uint64_t limit = 500;    // number of kmers to keep.
    float param_lc = 1.0;    // low complexity filter threshold, allow all known adapters to pass.
    uint64_t v = 1;          // verbosity
    bool skip_end = false;   // skip end adapter ressearch
    uint64_t nb_of_runs = 1; // Number of counts to run
    float lc = 1.0;          // adjusted LC filter threshold

    ////////////////////////
    // Config File 
    ////////////////////////

    getOptionValue(config_file, parser, "conf");
    // reading config file, if any, and adjusting variables.
    if(not config_file.empty() ){
        arg_map params = parse_config(config_file);
        param_lc  = params.count("lc" )>0 ? std::stof(params["lc"] ) : lc;
        k         = params.count("k"  )>0 ? std::stoi(params["k"]  ) : k;
        v         = params.count("v"  )>0 ? std::stoi(params["v"]  ) : v;
        sn        = params.count("sn" )>0 ? std::stoi(params["sn"] ) : sn;
        sl        = params.count("sl" )>0 ? std::stoi(params["sl"] ) : sl;
        limit     = params.count("lim")>0 ? std::stoi(params["lim"]) : limit;
        nb_thread = params.count("nt" )>0 ? std::stoi(params["nt"] ) : nb_thread;
        solid_km  = params.count("sk" )>0 ? std::stoi(params["sk"] ) : solid_km;
        skip_end  = params.count("se" )>0 ? true : false;
        forbid_kmer = params.count("fk") >0 ? params["fk"] : forbid_kmer;
        exact_out   = params.count("e")  >0 ? params["e"]  : exact_out;
        nb_of_runs  = params.count("mr" )>0 ? std::stoi(params["mr"] ) : nb_of_runs;    
    }

    /////////////////////////////////
    // Assigning values to variables.
    /////////////////////////////////

    // If options have been manually set, override config.
    getOptionValue(limit, parser, "lim");
    getOptionValue(param_lc, parser, "lc");
    getOptionValue(k, parser, "k");
    getOptionValue(v, parser, "v");
    getOptionValue(sl, parser, "sl");
    getOptionValue(sn, parser, "sn");
    getOptionValue(nb_thread, parser, "nt");
    getOptionValue(output, parser, "o");
    getOptionValue(exact_out, parser, "e");
    getOptionValue(forbid_kmer, parser, "fk");
    getOptionValue(solid_km, parser, "sk");
    getOptionValue(nb_of_runs, parser, "mr");

    // Except for flags, check if they are set in either config or manually
    skip_end = skip_end or isSet(parser, "skip_end");

    // Input file, always required
    std::string input_file;
    getArgumentValue(input_file, parser, 0);

    // Set of forbidden k-mers.
    kmer_set_t kmer_set;
    if(not forbid_kmer.empty()){
        print("Parsing the fobidden kmer list");
        kmer_set = parse_kmer_list(forbid_kmer);
    }

    // Adjusting "multi run" verbosity to avoid flooding stdout.
    int mr_v = v;
    if(nb_of_runs > 1 and v < 2){
        mr_v =  0;
    }

    std::string warning = "/!\\ WARNING: ";
    std::string error_pref = "/!\\ ERROR: ";

    // checking value for k
    if( k<2 or k>32 ){
        throw std::invalid_argument(error_pref + "kmer size must be between 2 and 32 (included)");
    }
    // checking if  k is bigger than sampling length
    if( k > sl ){
        throw std::invalid_argument(error_pref + "kmer size must be smaller than the sampling length (k <= sl)");
    }

    // adjusting low complexity to kmer size, from a k16 base.
    lc = adjust_threshold(param_lc, 16, k );

    // print parameters
    if(v>0){
        std::cout << "Kmer size:             " << k         << std::endl;
        std::cout << "Sampled sequences:     " << sn        << std::endl;
        std::cout << "Sampling length        " << sl        << std::endl;
        std::cout << "LC filter threshold:   " << param_lc  << std::endl;
        std::cout << "Adjusted LC threshold: " << lc        << std::endl;
        std::cout << "Nb thread:             " << nb_thread << std::endl;
        if(solid_km != 0){
            std::cout << "Solid kmers:           " << solid_km << std::endl;
        }
        else{
            std::cout << "Number of kept kmer:   " << limit     << std::endl;
        }
        std::cout << "Number of runs:        " << nb_of_runs << std::endl;
        std::cout << "Verbosity level:       " << v          << std::endl;
    }

    // number of tab to display
    int tab_level = 0;

    // if we run counts more than one time
    if(v > 0 and nb_of_runs >1){
        std::cout << "\nA total of " << nb_of_runs <<" runs will be performed."<< std::endl;
    }

    // Parsing input fasta file.
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    if(v>0){
        print("Parsing FASTA file",tab_level);
    }
    SeqFileIn seqFileIn(toCString(input_file));
    readRecords(ids, seqs, seqFileIn);

    // Display number of sequences
    if(v>0){
        unsigned nb_seq = length(seqs);
        print("Number of sequences found: " + std::to_string(nb_seq) + ".", tab_level);
    }

    // MAIN LOOP: STARTING K-MER COUNTING
    std::string run_suffix = "";
    for(uint64_t current_run = 0; current_run < nb_of_runs; current_run++){
        // If we run more than once, add a suffix
        run_suffix = "_" + std::to_string(current_run);
        if(nb_of_runs > 1){
            if(v>0)
                std::cout << "Starting run number " << current_run + 1 << std::endl;
        }
        // Checking if we can sample the requested number of sequences, else return the whole set
        uint64_t sequence_set_size = length(seqs);
        if(sn > sequence_set_size){ 
            std::cerr << warning << "Sequence set too small for the requested sample size\n";
            std::cerr << warning << "The whole set will be used.\n" ;
            sn = sequence_set_size;
        }

        // general flag for file output
        bool success = true;

        // performing ressearch on both ends
        std::array<std::string, 2 > ends = {"start","end"};

        bool bottom = false; // checking if we search top adapter(start) or bottom adapter (end)
        tab_level += 1;
        for(std::string which_end: ends){

            if(v>0){
                print("Working on sequence " + which_end + ".", tab_level -1);
            }
            if(mr_v>0){
                // sample and cut sequences to required length
                print("Sampling",tab_level);
            }
            sequence_set_type sample = sampleSequences(seqs, sn, sl, bottom, mr_v);


            // counting k-mers on the sampled sequences
            if(mr_v>0){
                print("Exact k-mer count", tab_level);
            }
            counter count = count_kmers(sample, k, lc, kmer_set);
            
            //  DEBUG
            /*pair_vector all_kmer = get_most_frequent(count, sequence_set_size, k);
            exportCounter(all_kmer, k, "all_kmers" + run_suffix + ".txt");
            */

            // keeping most frequents kmers
            if(mr_v>0){
                print("Number of kmer found: " + std::to_string(count.size()), tab_level);
            }
            
            // Either keep the most frequent kmers or keep solid kmers.
            pair_vector first_n_vector;
            if(solid_km != 0){
                if(mr_v>0){
                    print("Keeping solid k-mer", tab_level);
                }
                first_n_vector = get_solid_kmers(count, solid_km);
            }
            else{
                if(mr_v>0){
                    print("Keeping most frequent k-mer", tab_level);
                }
                first_n_vector = get_most_frequent(count, limit, k);
            }
            
            if(mr_v>0){
                print("Number of kmer kept:  " + std::to_string(first_n_vector.size()), tab_level ) ;
            }
            

            // Exporting exact kmer count, if required
            if(not exact_out.empty() ){
                if(mr_v>0)
                    print("Exporting exact kmer count", tab_level);
                success = exportCounter(first_n_vector, k, exact_out + run_suffix + "." + which_end );
                if(!success){
                    std::cerr << error_pref + "Failed to export exact k-mer count" << std::endl ;
                    std::cerr << "Path: " << exact_out + run_suffix + "." + which_end << std::endl ;
                    return(1);
                }
            }

            // Counting with at most 2 errors
            if(mr_v>0){
                print("Approximate k-mer count",tab_level);
            }
            counter error_counter = errorCount(sample, first_n_vector, nb_thread, k, mr_v);
            pair_vector sorted_error_count = get_most_frequent(error_counter, limit, k);

            if(mr_v>0){
                print("Exporting approximate count",tab_level);
            }
            success = exportCounter(sorted_error_count,k, output + run_suffix + "." + which_end );
            if(!success){
                    std::cerr << error_pref + "Failed to export approximate k-mer count" << std::endl ;
                    std::cerr << "Path: " << output + run_suffix + "." + which_end << std::endl ;
                    return(1);
            }


            if(mr_v>0){
                print("Done",tab_level);
            }
            
            clear(sample);
            
            // Shall we process read end ?
            if(skip_end){
                if(mr_v>0){
                    print("Skipping end adapter ressearch");
                    break;
                }
            }
            else{
                bottom = true;
            }
          
        } 
        tab_level --;
    }

    return 0;
}
