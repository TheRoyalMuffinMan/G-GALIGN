#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

#define DEFAULT_OUTPUT_FILE "output.txt"
#define DEFAULT_GPU_ID 0
#define DEFAULT_THREAD_SIZE 256

typedef struct {
    int64_t score;
    std::string updated_query;
    std::string alignment;
    std::string updated_ref;
} Result;

int write_results(int score, std::string filename, 
                  std::string reference_id, std::string reference, 
                  std::string alignment, std::string query_id, std::string query);
#endif