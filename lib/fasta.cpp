#include <iostream>
#include <iterator>
#include <fstream>
#include "fasta.hpp"

Fasta::Fasta(std::string path) {
    this->path = path;
}

Fasta::~Fasta() {

}

void Fasta::print() const {
    std::cout << "FASTA:" << std::endl;
    std::cout << "File: " << this->path << std::endl;
    for (size_t i = 0; i < genes.size(); i++) {
        std::cout << genes[i].id << ": " << genes[i].sequence << std::endl;
    }
    std::cout << std::endl;
}

int Fasta::read() {
    std::ifstream file(this->path);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Couldn't open or create file" << std::endl;
        return 1;
    }

    while (std::getline(file, line).good()) {
        // Empty line, ignore and continue
        if (line.empty()) {
            continue;
        }

        // Take 2 lines and store in unordered_map [metadata] => sequence
        if (line[0] == '>') {
            // Drop '>' character
            std::string metadata = line.substr(1);
            std::string sequence;

            // Read in sequence
            if (!std::getline(file, sequence)) {
                std::cerr << "Error: EOF before matching current metadata" << std::endl;
                return 1;
            }
            this->genes.push_back({metadata, sequence});
        } else {
            std::cerr << "Error: Doesn't match FASTA file format" << std::endl;
            return 1;
        }
    }

    file.close();
    return 0;
}