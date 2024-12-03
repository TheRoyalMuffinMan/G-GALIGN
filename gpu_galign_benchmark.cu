#include <iostream>
#include <cuda_runtime_api.h>
#include "include/cmd.hpp"
#include "include/fasta.hpp"
#include "include/globals.hpp"

// Set to 0 to disable debugging
#define DDEBUG 0

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


void print_matrix(int64_t *matrix, int64_t m, int64_t n) {
    int count = 0;
    std::cout << "DP Matrix: " << std::endl;
    for (int64_t i = 0; i < m; ++i) {
        for (int64_t j = 0; j < n; ++j) {
            std::cout << matrix[i * n + j] << " ";
            if (matrix[i * n + j] != 0) count++;
        }
        std::cout << std::endl;
    }
    std::cout << "Non-zero entries: " << count << std::endl;
}

__global__ void init_gaps(int64_t* dp_table, int64_t m, int64_t n, int gap_penalty) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Initializes the first column for all rows (dp_table[idx * n])
    if (idx < m) {
        int64_t row_position = idx * n;
        dp_table[row_position] = gap_penalty * idx;
        #if (DDEBUG != 0)
            printf("Initializing first column at position %lld (row %d, col 0)\n", row_position, idx);
        #endif
    }

    // Initializes the first row for all columns (dp_table[idx])
    if (idx < n) {
        dp_table[idx] = gap_penalty * idx;
        #if (DDEBUG != 0)
            printf("Initializing first row at position %d (row 0, col %d)\n", idx, idx);
        #endif
    }
}


Result needleman_wunsch(std::string reference, std::string query, 
                        int gap_penalty, int mismatch_penalty, 
                        int match_score, int ignore_outer_gaps) {
    Result res = {0};
    int64_t m = query.size() + 1, n = reference.size() + 1;
    size_t data_size = sizeof(int64_t) * m * n;
    int device_count;
    int64_t *device_dp_table = nullptr;
    int64_t *host_dp_table = new int64_t[m * n];

    // Get GPU counts
    checkCudaErrors(cudaGetDeviceCount(&device_count));

    // Utilize single GPU for processing
    checkCudaErrors(cudaSetDevice(DEFAULT_GPU_ID));

    // Allocate memory on the GPU for the DP matrix (might need to tile this)
    checkCudaErrors(cudaMalloc(&device_dp_table, data_size));

    // Initialize the gaps of the row and column
    int64_t nthreads = 256;
    int64_t nblocks = (std::max(m, n) + nthreads - 1) / nthreads;
    #if (DDEBUG != 0)
        std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
        std::cout << "(m = " << m << ", n = " << n << ")" << std::endl;
    #endif
    init_gaps<<<nthreads, nblocks>>>(device_dp_table, m, n, gap_penalty);

    // Begin needleman wunsch anti-diagonal approach on the GPU


    // Copy the result matrix from device to host
    checkCudaErrors(cudaMemcpy(host_dp_table, device_dp_table, data_size, cudaMemcpyDeviceToHost));

    // Print the 2D matrix (from the host memory)
    #if (DDEBUG != 0)
        print_matrix(host_dp_table, m, n);
    #endif

    // Free the device memory after use
    cudaFree(device_dp_table);
    delete[] host_dp_table;

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

    Result res = needleman_wunsch(
        reference.sequence,
        query.sequence,
        args.gap_penalty,
        args.mismatch_penalty,
        args.match_score,
        args.ignore_outer_gaps
    );

    write_results(10, args.output, reference.id, reference.sequence, "", query.id, query.sequence);
}