#ifndef FASTA_H
#define FASTA_H

#include <string>
#include <unordered_map>

class Fasta {
    public:
        Fasta(std::string path);
        ~Fasta();
        void print() const;
        int read();

        std::unordered_map<std::string, std::string> genes;
    private:
        std::string path;
};

#endif