#include <iostream>
#include <fstream>
#include "globals.hpp"

int write_results(int score, std::string filename, 
                 std::string reference_id, std::string reference, 
                 std::string alignment, std::string query_id, std::string query) {               
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Couldn't open or create file" << std::endl;
        return 1;
    }

    file << "Alignment Score: " << score << std::endl;
    file << ">" << reference_id << std::endl;
    file << reference << std::endl;
    file << alignment << std::endl;
    file << query << std::endl;
    file << ">" << query_id << std::endl;

    file.close();
    return 0;
}