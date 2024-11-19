#include <iostream>
#include "include/cmd.hpp"
#include "include/parser.hpp"

// Comment out to disable debugging
#define ENABLE_DEBUG 

int main(int argc, char *argv[]) {
    CommandLineArgs args(argc, argv);
    if (args.parse() != 0) {
        std::exit(EXIT_FAILURE);
    }

    #ifdef ENABLE_DEBUG
        args.print();
    #endif

}