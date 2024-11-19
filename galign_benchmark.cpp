#include <iostream>
#include "include/cmd.hpp"
#include "include/parser.hpp"

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

}