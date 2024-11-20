#include <iostream>
#include "include/cmd.hpp"
#include "include/fasta.hpp"

// Set to 0 to disable debugging
#define DDEBUG 1

int main(int argc, char *argv[]) {
    CommandLineArgs args(argc, argv);
    if (args.parse() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        args.print();
    #endif

    Fasta query(args.query);
    Fasta reference(args.reference);
    if (query.read() != 0 || reference.read() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        query.print();
        reference.print();
    #endif

}