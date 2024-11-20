#include <iostream>
#include <getopt.h>
#include "globals.hpp"
#include "cmd.hpp"

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

void CommandLineArgs::print() const {
    std::cout << "Command Line Arguments:" << std::endl;
    std::cout << "Query: " << this->query << std::endl;
    std::cout << "Reference: " << this->reference << std::endl;
    std::cout << "Output: " << this->output << std::endl;
    std::cout << "Gap Penalty: " << this->gap_penalty << std::endl;
    std::cout << "Mismatch Penalty: " << this->mismatch_penalty << std::endl;
    std::cout << "Match Score: " << this->match_score << std::endl;
    std::cout << "Ignore Outer Gaps: " << (this->ignore_outer_gaps ? "true" : "false") << std::endl;
    std::cout << std::endl;
}

int CommandLineArgs::parse() {
    // Setup flags, struct for arguments, and argument flag character
    int query_flag = 0, reference_flag = 0, output_flag = 0;
    int gap_flag = 0, mismatch_flag = 0, match_flag = 0;
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

    // Set defaults for optional arguments
    this->ignore_outer_gaps = 0;

    while ((ch = getopt_long(this->argc, this->argv, "q:r:o:g:p:m:i:h", long_options, NULL)) != -1) {
        switch (ch) {
            case 'q':
                this->query = std::string(optarg);
                query_flag = 1;
                break;
            case 'r':
                this->reference = std::string(optarg);
                reference_flag = 1;
                break;
            case 'o':
                this->output = std::string(optarg);
                output_flag = 1;
                break;
            case 'g':
                this->gap_penalty = std::stoi(optarg);
                gap_flag = 1;
                break;
            case 'p':
                this->mismatch_penalty = std::stoi(optarg);
                mismatch_flag = 1;
                break;
            case 'm':
                this->match_score = std::stoi(optarg);
                match_flag = 1;
                break;
            case 'i':
                this->ignore_outer_gaps = !!std::stoi(optarg);
                break;
            case 'h':
                goto error;
            default:
                std::cout << "Invalid option or missing argument" << std::endl;
                goto error;
        }
    }

    if (query_flag == 0) {
        std::cerr << "Error: No query file passed in" << std::endl;
        goto error;
    }

    if (reference_flag == 0) {
        std::cerr << "Error: No reference file passed in" << std::endl;
        goto error;
    }

    if (output_flag == 0) {
        std::cerr << "Warning: No output file passed in... defaulting to output.txt as filename" << std::endl;
        this->output = DEFAULT_OUTPUT_FILE;
    }

    if (gap_flag == 0) {
        std::cerr << "Error: No gap penalty passed in" << std::endl;
        goto error;
    }

    if (mismatch_flag == 0) {
        std::cerr << "Error: No mismatch penalty passed in" << std::endl;
        goto error;
    }

    if (match_flag == 0) {
        std::cerr << "Error: No match score passed in" << std::endl;
        goto error;
    }

    return 0;
error:
    usage();
    return 1;
}