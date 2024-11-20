#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

#define DEFAULT_OUTPUT_FILE "output.txt"

int write_results(int score, std::string filename, 
                  std::string reference_id, std::string reference, 
                  std::string alignment, std::string query_id, std::string query);
#endif