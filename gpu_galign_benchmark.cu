#include <iostream>
#include <chrono>
#include <algorithm>
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

// Prints the matrix (USED FOR DEBUG PURPOSES)
void print_matrix(int64_t *table, int64_t m, int64_t n) {
    std::cout << "DP Matrix: " << std::endl;
    for (int64_t i = 0; i < m; ++i) {
        for (int64_t j = 0; j < n; ++j) {
            std::cout << table[i * n + j] << " ";
        }
    }
}

// Helps compute the 1-D position using 2-D indexing for the table
__host__ __device__ int get_position(int64_t row, int64_t col, int64_t num_cols) {
    return row * num_cols + col;
}

// Initializes the gap penalties for the first row and column
__global__ void init_gaps(int64_t *dp_table, int64_t m, int64_t n, int gap_penalty) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    // Initializes the first column for all rows (dp_table[idx * n])
    if (idx < m) {
        int64_t row_position = idx * n;
        dp_table[row_position] = gap_penalty * idx;
    }

    // Initializes the first row for all columns (dp_table[idx])
    if (idx < n) {
        dp_table[idx] = gap_penalty * idx;
    }
}

// Processes an anti-diagonal in the DP table using start and end positions
__global__ void solve_anti_diagonal(char *query, char *reference,
                                    int64_t *dp_table, size_t data_size, 
                                    int64_t m, int64_t n, 
                                    int64_t start_row, int64_t start_col,
                                    int64_t end_row, int64_t end_col,
                                    int gap_penalty, int mismatch_penalty, int match_score) {
                                        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int64_t new_row = start_row - idx, new_col = start_col + idx;
    int64_t position = get_position(new_row, new_col, n);

    // Checks if the thread is within range and we aren't at the end position
    if (position < data_size && new_row >= end_row && new_col < end_col) {
        int64_t diagonal = 0, top = 0, left = 0;
        if (query[new_row - 1] == reference[new_col - 1]) {
            diagonal = dp_table[get_position(new_row - 1, new_col - 1, n)] + match_score;
        } else {
            diagonal = dp_table[get_position(new_row - 1, new_col - 1, n)] + mismatch_penalty;
        }

        top = dp_table[get_position(new_row - 1, new_col, n)] + gap_penalty;
        left = dp_table[get_position(new_row, new_col - 1, n)] + gap_penalty;

        dp_table[position] = max(diagonal, max(top, left));
    }
}

// Generates the traceback strings and acquires alignment score for the CPU
// NOTE: This is done serially (or with 1 GPU thread)
__global__ void build_traceback(char *query, char *reference, 
                                int64_t *dp_table, int64_t m, int64_t n, int64_t *alignment_score,
                                char *updated_query, char *alignment, char *updated_reference,
                                int64_t match_score, int64_t mismatch_penalty, int64_t gap_penalty) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx != 0) return; // Single-threaded traceback

    // Save alignment score for these two sequences
    *alignment_score = dp_table[get_position(m - 1, n - 1, n)];

    // Build mapping from computed DP table
    int64_t row = m - 1, col = n - 1;
    int64_t q_i = 0, a_i = 0, r_i = 0;
    while (row > 0 && col > 0) {
        int64_t match = dp_table[get_position(row - 1, col - 1, n)] + match_score;
        int64_t mismatch = dp_table[get_position(row - 1, col - 1, n)] + mismatch_penalty;
        int64_t top = dp_table[get_position(row - 1, col, n)] + gap_penalty;
        int64_t left = dp_table[get_position(row, col - 1, n)] + gap_penalty;

        // Two sequences matched in this position
        if (query[row - 1] == reference[col - 1] && match == dp_table[get_position(row, col, n)]) {
            alignment[a_i++] = '|';
            updated_query[q_i++] = query[row - 1];
            updated_reference[r_i++] = reference[col - 1];
            row = row - 1;
            col = col - 1;
        } 
        // Two sequences mismatched in this position
        else if (mismatch == dp_table[get_position(row, col, n)]) {
            alignment[a_i++] = 'x';
            updated_query[q_i++] = query[row - 1];
            updated_reference[r_i++] = reference[col - 1];
            row = row - 1;
            col = col - 1;
        } 
        // Came from top gap to get to this position
        else if (top == dp_table[get_position(row, col, n)]) {
            alignment[a_i++] = ' ';
            updated_query[q_i++] = query[row - 1];
            updated_reference[r_i++] = '_';
            row = row - 1;
        } 
        // Came from left gap to get to this position
        else if (left == dp_table[get_position(row, col, n)]) {
            alignment[a_i++] = ' ';
            updated_query[q_i++] = '_';
            updated_reference[r_i++] = reference[col - 1];
            col = col - 1;
        }
    }

    // Clean up remaining rows, move cells upward until we reach (0, 0)
    while (row > 0) {
        alignment[a_i++] = ' ';
        updated_query[q_i++] = query[row - 1];
        updated_reference[r_i++] = '_';
        row = row - 1;
    }

    // Clean up remaining columns, move cells leftward until we reach (0, 0)
    while (col > 0) {
        alignment[a_i++] = ' ';
        updated_query[q_i++] = '_';
        updated_reference[r_i++] = reference[col - 1];
        col = col - 1;
    }
    
    // Set null terminating character to end aligned strings
    alignment[a_i] = '\0';
    updated_query[q_i] = '\0';
    updated_reference[r_i] = '\0';
}

// Computes the needleman wunsch sequence alignment algorithm on the GPU
Result parallel_needleman_wunsch(std::string reference, std::string query, 
                                 int gap_penalty, int mismatch_penalty, int match_score) {

    Result res = {.score = 0, .updated_query = "", .alignment = "", .updated_ref = ""};
    int64_t m = query.size() + 1, n = reference.size() + 1;
    int64_t anti_diagonals = m + n - 3; // Discount 2 gap (row and column)
    int64_t nthreads;
    int64_t nblocks;
    size_t data_size = sizeof(int64_t) * m * n;
    size_t alignment_size = m + n + 1;
    int device_count;

    // Device (GPU) Addresses
    int64_t *device_dp_table = nullptr;
    int64_t *device_alignment_score;
    char *device_query = nullptr, *device_reference = nullptr;
    char *device_updated_query = nullptr, *device_alignment = nullptr, *device_updated_reference = nullptr;

    // Host (CPU) Addresses
    int64_t *host_alignment_score;
    char *host_updated_query = nullptr, *host_alignment = nullptr, *host_updated_reference = nullptr;

    // Begin timing for memory operations + DP table computation
    auto memory_execution_start = std::chrono::high_resolution_clock::now();

    // Get GPU counts
    checkCudaErrors(cudaGetDeviceCount(&device_count));

    // Utilize single GPU for processing
    checkCudaErrors(cudaSetDevice(DEFAULT_GPU_ID));

    // Allocate pinned (page-locked) memory on the host
    checkCudaErrors(cudaHostAlloc(&host_alignment_score, sizeof(int64_t), cudaHostAllocDefault));
    checkCudaErrors(cudaHostAlloc(&host_updated_query, alignment_size, cudaHostAllocDefault));
    checkCudaErrors(cudaHostAlloc(&host_alignment, alignment_size, cudaHostAllocDefault));
    checkCudaErrors(cudaHostAlloc(&host_updated_reference, alignment_size, cudaHostAllocDefault));

    // Allocate memory on the GPU for the DP table
    checkCudaErrors(cudaMalloc(&device_dp_table, data_size));

    // Allocate memory on the GPU for alignment score
    checkCudaErrors(cudaMalloc(&device_alignment_score, sizeof(int64_t)));

    // Allocate memory for query and reference strings on the GPU
    checkCudaErrors(cudaMalloc(&device_query, query.size()));
    checkCudaErrors(cudaMalloc(&device_reference, reference.size()));

    // Over-allocate memory for the aligned query, alignment, and aligned reference string 
    // (the strings shouldn't be greater than m + n)
    checkCudaErrors(cudaMalloc(&device_updated_query, alignment_size));
    checkCudaErrors(cudaMalloc(&device_alignment, alignment_size));
    checkCudaErrors(cudaMalloc(&device_updated_reference, alignment_size));

    // Copy query and reference strings to the GPU
    checkCudaErrors(cudaMemcpy(device_query, query.c_str(), query.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(device_reference, reference.c_str(), reference.size(), cudaMemcpyHostToDevice));

    // Get start time for computing table (include initialize gaps and compute DP)
    auto execution_start = std::chrono::high_resolution_clock::now();

    // Initialize the gaps of the row and column
    nthreads = DEFAULT_THREAD_SIZE;
    nblocks = (std::max(m, n) + nthreads - 1) / nthreads;
    #if (DDEBUG != 0)
        std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
        std::cout << "(m = " << m << ", n = " << n << ")" << std::endl;
    #endif
    init_gaps<<<nthreads, nblocks>>>(device_dp_table, m, n, gap_penalty);

    /****   Begin needleman wunsch anti-diagonal approach on the GPU   ****/
    /***  We will spawn out a kernel invocation for each anti-diagonal  ***/
    int64_t start_row = 0, start_col = 1;
    int64_t end_row = 0, end_col = 1;
    for (int64_t diag = 0; diag < anti_diagonals; diag++) {
        
        // Computes the start position and end position of the given
        // anti-diagonal in order to process it on the GPU, ex:
        /*
         *         e   e   e   e
         *         V   V   V   V
         * 0   0   0   0   0   0
         * 6 s>1   2   3   4   5 < e
         * 0 s>1   2   3   4   5 < e
         * 0 s>1   2   3   4   5 < e
         * 0 s>1   2   3   4   5 < e
         * 0 s>1   2   3   4   5 < e
         *         Ʌ   Ʌ   Ʌ   Ʌ
         *         s   s   s   s
        */
        int64_t cells = 0;
        if (start_row < m - 1) {
            start_row++;
            if (end_col + 1 < n) {
                end_col++, cells = start_row - end_row;
            } else {
                end_col = n, end_row++, cells = end_col - start_col;
            }
        } else {
            start_col++;
            if (end_col + 1 < n) {
                end_col++, cells = start_row - end_row;
            } else {
                end_col = n, end_row++, cells = end_col - start_col;
            }
        }

        // Spawn out threads based on distance and solve current anti-diagonal
        // Uses the number of cells (distance between start and end positions) to determine number of threads and blocks
        nthreads = std::min(static_cast<int64_t>(DEFAULT_THREAD_SIZE), cells);
        nblocks = (cells + nthreads - 1) / nthreads;
        #if (DDEBUG != 0)
            std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
            std::cout << "start_row, start_col: (" << start_row << ", " << start_col 
                      << "), end_row, end_col: (" << end_row << ", " << end_col << ")" 
                      << std::endl;
        #endif
        solve_anti_diagonal<<<nthreads, nblocks>>>(
            device_query, device_reference,
            device_dp_table, data_size,
            m, n, 
            start_row, start_col, 
            end_row, end_col,
            gap_penalty, mismatch_penalty, match_score
        );
    }

    // Get end time for computing table (include initialize gaps and compute DP)
    auto execution_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> execution_duration = execution_end - execution_start;
    std::cout << "Execution time: " 
              << execution_duration.count() << " seconds" << std::endl;

    // Builds the traceback computed from needleman wunsch using single thread and block
    build_traceback<<<1, 1>>>(
        device_query, device_reference,
        device_dp_table, m, n, device_alignment_score, 
        device_updated_query, device_alignment, device_updated_reference,
        match_score, mismatch_penalty, gap_penalty
    );

    // Copy score, updated query, alignment, and updated reference strings from the GPU back to the CPU
    checkCudaErrors(cudaMemcpy(host_alignment_score, device_alignment_score, sizeof(int64_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(host_updated_query, device_updated_query, alignment_size, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(host_alignment, device_alignment, alignment_size, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(host_updated_reference, device_updated_reference, alignment_size, cudaMemcpyDeviceToHost));
    
    // End timing for memory operations + DP table computation
    auto memory_execution_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> memory_execution_duration = memory_execution_end - memory_execution_start;
    std::cout << "Memory + Execution time: " 
              << memory_execution_duration.count() << " seconds" << std::endl;
    
    // Get score and convert character arrays to strings
    res.score = *host_alignment_score;
    res.updated_query = std::string(host_updated_query);
    res.alignment = std::string(host_alignment);
    res.updated_ref = std::string(host_updated_reference);

    // Free the host memory after use
    cudaFreeHost(host_updated_query);
    cudaFreeHost(host_alignment);
    cudaFreeHost(host_updated_reference);
    cudaFreeHost(host_alignment_score);

    // Free the device memory after use
    cudaFree(device_query);
    cudaFree(device_reference);
    cudaFree(device_dp_table);
    cudaFree(device_updated_query);
    cudaFree(device_alignment);
    cudaFree(device_updated_reference);
    cudaFree(device_alignment_score);

    // Reverse the sequences for the output file
    std::reverse(res.updated_query.begin(), res.updated_query.end());
    std::reverse(res.alignment.begin(), res.alignment.end());
    std::reverse(res.updated_ref.begin(), res.updated_ref.end());

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

    Result res = parallel_needleman_wunsch(
        reference.sequence,
        query.sequence,
        args.gap_penalty,
        args.mismatch_penalty,
        args.match_score
    );

    write_results(
        res.score, args.output, 
        reference.id, res.updated_ref, 
        res.alignment, 
        query.id, res.updated_query
    );
}