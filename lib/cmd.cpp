#include "cmd.h"
#include <iostream>
#include <getopt.h>

CommandLineArgs::CommandLineArgs(int argc, char* argv[]) {
    this->argc = argc;
    this->argv = argv;
}

CommandLineArgs::~CommandLineArgs() {}

void CommandLineArgs::usage() const {
    std::cout << "Usage: program_name [OPTIONS]\n"
              << "Options:\n"
              << "  -q, --query <file>              Input query file\n"
              << "  -r, --reference <file>          Input reference file\n"
              << "  -o, --output <file>             Output file\n"
              << "  -g, --gap_penalty <int>         Gap penalty\n"
              << "  -p, --mismatch_penalty <int>    Mismatch penalty\n"
              << "  -m, --match_score <int>         Match score\n"
              << "  -i, --ignore_outer_gaps         Ignore outer gaps\n"
              << "  -h, --help                      Show this help message and exit\n";
}

void CommandLineArgs::parse() {
    struct option long_options[] = {
        {"query", required_argument, NULL, 'q'},
        {"reference", required_argument, NULL, 'r'},
        {"output", required_argument, NULL, 'o'},
        {"gap_penalty", required_argument, NULL, 'g'},
        {"mismatch_penalty", required_argument, NULL, 'p'},
        {"match_score", required_argument, NULL, 'm'},
        {"ignore_outer_gaps", optional_argument, 0, 'i'},
        {"help", no_argument, NULL, 'h'}
    };

    char ch;
    while ((ch = getopt_long(this->argc, this->argv, "q:r:o:g:p:m:i:h", long_options, NULL))) {
        switch (ch) {
            case 'q':
                this->query = std::string(optarg);
                break;
            case 'r':
                this->reference = std::string(optarg);
                break;
            case 'o':
                this->output = std::string(optarg);
                break;
            case 'g':
                this->gap_penalty = std::stoi(optarg);
                break;
            case 'p':
                this->mismatch_penalty = std::stoi(optarg);
                break;
            case 'm':
                this->match_score = std::stoi(optarg);
                break;
            case 'i':
                this->ignore_outer_gaps = !!std::stoi(optarg);
                break;
            case 'h':
                usage();
            default:
                std::cout << "Invalid option or missing argument" << std::endl;
                usage();
        }
    }
}
