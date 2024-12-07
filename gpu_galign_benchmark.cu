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

// Helps compute the 1-D position using 2-D indexing for the matrix
__host__ __device__ int get_position(int64_t row, int64_t col, int64_t num_cols) {
    return row * num_cols + col;
}

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

__global__ void solve_anti_diagonals(char *query, char *reference,
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

Result needleman_wunsch(std::string reference, std::string query, 
                        int gap_penalty, int mismatch_penalty, 
                        int match_score, int ignore_outer_gaps) {
    Result res = {.score = 0, .updated_query = "", .alignment = "", .updated_ref = ""};
    int64_t m = query.size() + 1, n = reference.size() + 1;
    int64_t *device_dp_table = nullptr;
    int64_t *host_dp_table = new int64_t[m * n];
    int64_t anti_diagonals = m + n - 3; // Discount 2 gap (row and column)
    size_t data_size = sizeof(int64_t) * m * n;
    char *device_query = nullptr, *device_reference = nullptr;
    int device_count;

    // Get GPU counts
    checkCudaErrors(cudaGetDeviceCount(&device_count));

    // Utilize single GPU for processing
    checkCudaErrors(cudaSetDevice(DEFAULT_GPU_ID));

    // Allocate memory on the GPU for the DP matrix (might need to tile this)
    checkCudaErrors(cudaMalloc(&device_dp_table, data_size));

    // Allocate memory for query and reference strings on the GPU
    checkCudaErrors(cudaMalloc(&device_query, query.size()));
    checkCudaErrors(cudaMalloc(&device_reference, reference.size()));

    // Copy query and reference strings to GPU
    checkCudaErrors(cudaMemcpy(device_query, query.c_str(), query.size(), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(device_reference, reference.c_str(), reference.size(), cudaMemcpyHostToDevice));

    // Initialize the gaps of the row and column
    auto start = std::chrono::high_resolution_clock::now();
    int64_t nthreads = DEFAULT_THREAD_SIZE;
    int64_t nblocks = (std::max(m, n) + nthreads - 1) / nthreads;
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
         *         ^   ^   ^   ^
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
        // Uses the number of cells (distance between start and end positions) to determine
        // number of threads and blocks
        int64_t nthreads = std::min(static_cast<int64_t>(DEFAULT_THREAD_SIZE), cells);
        int64_t nblocks = (cells + nthreads - 1) / nthreads;
        solve_anti_diagonals<<<nthreads, nblocks>>>(
            device_query, device_reference,
            device_dp_table, data_size,
            m, n, 
            start_row, start_col, 
            end_row, end_col,
            gap_penalty, mismatch_penalty, match_score
        );
    }

    // Insure all kernels have finished executing
    cudaDeviceSynchronize();

    // Get time taken to compute matrix (include initialize gaps and compute DP)
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " 
              << duration.count() << " seconds" << std::endl;  // Convert seconds to microseconds
    
    // Copy the result matrix from device to host
    checkCudaErrors(cudaMemcpy(host_dp_table, device_dp_table, data_size, cudaMemcpyDeviceToHost));

    // Print the 2D matrix (from the host memory)
    #if (DDEBUG != 0)
        print_matrix(host_dp_table, m, n);
    #endif

    // Build mapping from computed DP table from the GPU
    int64_t row = m - 1, col = n - 1;
    res.score = host_dp_table[get_position(row, col, n)];
    while (row > 0 && col > 0) {
        int64_t match = host_dp_table[get_position(row - 1, col - 1, n)] + match_score;
        int64_t mismatch = host_dp_table[get_position(row - 1, col - 1, n)] + mismatch_penalty;
        int64_t top = host_dp_table[get_position(row - 1, col, n)] + gap_penalty;
        int64_t left = host_dp_table[get_position(row, col - 1, n)] + gap_penalty;

        // Two sequences matched in this position
        if (query[row - 1] == reference[col - 1] && match == host_dp_table[get_position(row, col, n)]) {
            res.alignment += "|";
            res.updated_query += query[row - 1];
            res.updated_ref += reference[col - 1];
            row = row - 1, col = col - 1;
        
        // Two sequences mismatched in this position
        } else if (mismatch == host_dp_table[get_position(row, col, n)]) {
            res.alignment += "x";
            res.updated_query += query[row - 1];
            res.updated_ref += reference[col - 1];
            row = row - 1, col = col - 1;
        
        // Came from top gap to get to this position
        } else if (top == host_dp_table[get_position(row, col, n)]) {
            res.alignment += " ";
            res.updated_query += query[row - 1];
            res.updated_ref += "_";
            row = row - 1;
        
        // Came from left gap to get to this position
        } else if (left == host_dp_table[get_position(row, col, n)]) {
            res.alignment += " ";
            res.updated_query += "_";
            res.updated_ref += reference[col - 1];
            col = col - 1;
        }
    }

    // Clean up remaining rows, move cells upward until we reach (0, 0)
    while (row > 0) {
        res.alignment += " ";
        res.updated_query += query[row - 1];
        res.updated_ref += "_";
        row = row - 1;
    }

    // Clean up remaining rows, move cells leftward until we reach (0, 0)
    while (col > 0) {
        res.alignment += " ";
        res.updated_query += "_";
        res.updated_ref += reference[col - 1];
        col = col - 1;
    }

    // Free the device memory after use
    cudaFree(device_query);
    cudaFree(device_reference);
    cudaFree(device_dp_table);

    // Free allocated host memory
    delete[] host_dp_table;

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

    Result res = needleman_wunsch(
        reference.sequence,
        query.sequence,
        args.gap_penalty,
        args.mismatch_penalty,
        args.match_score,
        args.ignore_outer_gaps
    );

    write_results(
        res.score, args.output, 
        reference.id, res.updated_ref, 
        res.alignment, 
        query.id, res.updated_query
    );
}