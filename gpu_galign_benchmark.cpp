#include <iostream>
#include "include/cmd.hpp"
#include "include/fasta.hpp"
#include "include/globals.hpp"

// Set to 0 to disable debugging
#define DDEBUG 1

// CUDA Helper Directives
#define checkCudaErrors(call)                                           \
    do {                                                                \
        cudaError_t err = call;                                         \
        if (err != cudaSuccess) {                                       \
            std::cerr << "CUDA error at " << __FILE__ << " "            \
                      << __LINE__ << ": " << cudaGetErrorString(err)    \
                      << std::endl;                                     \
            std::exit(EXIT_FAILURE);                                    \
        }                                                               \
    } while (0)


// 1,073,741,824
Result needleman_wunsch(std::string reference, std::string query, 
                        int gap_penalty, int mismatch_penalty, 
                        int match_score, int ignore_outer_gaps) {
    Result res = {0};
    int *deviceDPTable = nullptr;
    int m = query.size() + 1, n = reference.size() + 1;
    int deviceCount;
    int deviceID;

    // Get GPU counts
    checkCudaErrors(cudaGetDeviceCount(&devicesCount));

    // Utilize single GPU for processing
    checkCudaErrors(cudaSetDevice(DEFAULT_GPU_ID));

    // Allocate memory on the GPU for the DP matrix (might need to tile this)
    checkCudaErrors(cudaMalloc(&deviceDPTable, sizeof(int) * m * n));



    return res;
}

int main(int argc, char *argv[]) {
    CommandLineArgs args(argc, argv);
    if (args.parse() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        args.print();
    #endif

    Fasta query_fasta(args.query);
    Fasta reference_fasta(args.reference);
    if (query_fasta.read() != 0 || reference_fasta.read() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        query_fasta.print();
        reference_fasta.print();
    #endif

    // Only consider the first entry in each file for global alignment
    Gene query = query_fasta.genes[0];
    Gene reference = reference_fasta.genes[0];



    write_results(10, args.output, reference.id, reference.sequence, "", query.id, query.sequence);
}