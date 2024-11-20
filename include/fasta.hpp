#ifndef FASTA_H
#define FASTA_H

#include <string>
#include <vector>

typedef struct {
    std::string id;
    std::string sequence;
} Gene;

class Fasta {
    public:
        Fasta(std::string path);
        ~Fasta();
        void print() const;
        int read();
        
        std::vector<Gene> genes;
    private:
        std::string path;
};

#endif