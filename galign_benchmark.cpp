#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>

typedef struct arguments {
    std::string query;
    std::string reference;
    std::string output;
    int gap_penalty;
    int mismatch_penalty;
    int match_score;
    bool ignore_outer_gaps;
} arguments_t;

arguments_t parse_arguments(int argc, char *argv[]) {
    struct option long_options[] = {
        {"query", required_argument, NULL, 'q'},
        {"reference", required_argument, NULL, 'r'},
        {"output", required_argument, NULL, 'o'},
        {"gap_penalty", required_argument, NULL, 'g'},
        {"mismatch_penalty", required_argument, NULL, 'p'},
        {"match_score", required_argument, NULL, 'm'},
        {"ignore_outer_gaps", optional_argument, 0, 'i'}
    };
    arguments_t args;
    char ch;

    while ((ch = getopt_long(argc, argv, "q:r:o:g:p:m:i::", long_options, NULL))) {
        switch (ch) {
            case 'q':
                args.query = std::string(optarg);
                break;
            case 'r':
                args.reference = std::string(optarg);
                break;
            case 'o':
                args.output = std::string(optarg);
                break;
            case 'g':
                args.gap_penalty = std::stoi(optarg);
                break;
            case 'p':
                args.mismatch_penalty = std::stoi(optarg);
                break;
            case 'm':
                args.match_score = std::stoi(optarg);
                break;
            case 'i':
                args.ignore_outer_gaps = !!std::stoi(optarg);
        }
    }

    return args;
}

int main(int argc, char *argv[]) {
    arguments_t args = parse_arguments(argc, argv);
    
}