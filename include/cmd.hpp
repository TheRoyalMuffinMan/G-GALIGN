#ifndef CMD_H
#define CMD_H

#include <string>

class CommandLineArgs {
    public:
        CommandLineArgs(int argc, char* argv[]);
        ~CommandLineArgs();
        void usage() const;
        void print() const;
        int parse();
    
        std::string query;
        std::string reference;
        std::string output;
        int gap_penalty;
        int mismatch_penalty;
        int match_score;
    private:
        int argc;
        char** argv;
};

#endif