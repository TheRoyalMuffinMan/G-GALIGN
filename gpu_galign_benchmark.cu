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

__host__ __device__ int get_position(int64_t row, int64_t col, int64_t num_cols) {
    return row * num_cols + col;
}

__global__ void init_gaps(int64_t *dp_table, int64_t m, int64_t n, int gap_penalty) {
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

__global__ void solve_anti_diagonals(char *query, char *reference,
                                     int64_t *dp_table, size_t data_size, 
                                     int64_t m, int64_t n, 
                                     int64_t start_row, int64_t start_col,
                                     int64_t end_row, int64_t end_col,
                                     int gap_penalty, int mismatch_penalty, int match_score) {
                                        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int64_t new_row = start_row - idx, new_col = start_col + idx;
    int64_t position = get_position(new_row, new_col, n);
    #if (DDEBUG != 0)
        printf("start_row: %lld, start_col: %lld, idx: %d, new_row: %lld, new_col: %lld, position: %lld\n", start_row, start_col, idx, new_row, new_col, position);
    #endif

    if (position == 7647 || (new_row == 2 && new_col == 1)) {
        printf("m: %lld, n: %lld\n", m, n);
        printf("start_row: %lld, start_col: %lld, end_row: %lld, end_col: %lld, threadidx: %d\n", 
            start_row, start_col, end_row, end_col, idx);
        printf("new_row: %lld, new_col: %lld, position: %lld\n", new_row, new_col, position);    
    }

    // Valid index
    if (position < data_size && new_row > end_row && new_col < end_col) {
        int64_t diagonal = 0, top = 0, left = 0;
        if (query[new_row - 1] == reference[new_col - 1]) {
            diagonal = dp_table[get_position(new_row - 1, new_col - 1, n)] + match_score;
        } else {
            diagonal = dp_table[get_position(new_row - 1, new_col - 1, n)] + mismatch_penalty;
        }

        top = dp_table[get_position(new_row - 1, new_col, n)] + gap_penalty;
        left = dp_table[get_position(new_row, new_col - 1, n)] + gap_penalty;

        int64_t max_value = diagonal;  // Start with the diagonal value

        if (top > max_value) {
            max_value = top;  // Update max_value if top is greater
        }

        if (left > max_value) {
            max_value = left;  // Update max_value if left is greater
        }

        // Finally assign the maximum value to dp_table[position]
        dp_table[position] = max_value;
    }
}

Result needleman_wunsch(std::string reference, std::string query, 
                        int gap_penalty, int mismatch_penalty, 
                        int match_score, int ignore_outer_gaps) {
    Result res = {.score = 0, .updated_query = "", .alignment = "", .updated_ref = ""};
    int64_t m = query.size() + 1, n = reference.size() + 1;
    size_t data_size = sizeof(int64_t) * m * n;
    int64_t *device_dp_table = nullptr;
    int64_t *host_dp_table = new int64_t[m * n];
    char *device_query = nullptr, *device_reference = nullptr;
    int device_count;

    // auto start = std::chrono::high_resolution_clock::now();
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
    int64_t nthreads = 256;
    int64_t nblocks = (m + n + nthreads - 1) / nthreads;
    #if (DDEBUG != 0)
        std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
        std::cout << "(m = " << m << ", n = " << n << ")" << std::endl;
    #endif
    auto start = std::chrono::high_resolution_clock::now();
    init_gaps<<<nthreads, nblocks>>>(device_dp_table, m, n, gap_penalty);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " 
              << duration.count() << " seconds" << std::endl;  // Convert seconds to microseconds

    /** Begin needleman wunsch anti-diagonal approach on the GPU **/
    /* We will spawn out a kernel invocation for each anti-diagonal */

    // Compute the first half
    int64_t start_col = 1, end_row = 0;
    int64_t end_col = 2;
    for (int64_t start_row = 1; start_row < m; start_row++) {
        // Initialization
        int64_t distance = std::min(start_row, end_col);

        // Setting up number of threads and block size
        int64_t nthreads = std::min(static_cast<int64_t>(256), distance);
        int64_t nblocks = (distance + nthreads - 1) / nthreads;
        #if (DDEBUG != 0)
            std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
            std::cout << "(m = " << m << ", n = " << n << ")" << std::endl;
        #endif
        solve_anti_diagonals<<<nthreads, nblocks>>>(
            device_query, device_reference,
            device_dp_table, data_size,
            m, n, 
            start_row, start_col, 
            end_row, end_col,
            gap_penalty, mismatch_penalty, match_score
        );
        end_col = std::min(end_col + 1, n);
    }

    // Compute the second half
    int64_t start_row = m - 1;
    end_col = n;
    for (start_col = 2; start_col < n; start_col++) {
        // Initialization
        end_row = start_col;
        int64_t distance = n - start_col;

        // Setting up number of threads and block size
        int64_t nthreads = std::min(static_cast<int64_t>(256), distance);
        int64_t nblocks = (distance + nthreads - 1) / nthreads;
        #if (DDEBUG != 0)
            std::cout << "(nThreads = " << nthreads << ", nBlocks = " << nblocks << ")" << std::endl;
            std::cout << "(m = " << m << ", n = " << n << ")" << std::endl;
        #endif
        solve_anti_diagonals<<<nthreads, nblocks>>>(
            device_query, device_reference,
            device_dp_table, data_size,
            m, n, 
            start_row, start_col, 
            end_row, end_col,
            gap_penalty, mismatch_penalty, match_score
        );
    }

    // Copy the result matrix from device to host
    checkCudaErrors(cudaMemcpy(host_dp_table, device_dp_table, data_size, cudaMemcpyDeviceToHost));

    // Get time taken to compute matrix
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // std::cout << "Execution time: " 
    //           << duration.count() << " seconds" << std::endl;  // Convert seconds to microseconds

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