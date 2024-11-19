#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include "include/cmd.hpp"

int main(int argc, char *argv[]) {
    CommandLineArgs args(argc, argv);
    args.parse();
}