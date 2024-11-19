#ifndef CMD_H
#define CMD_H

#include <iostream>
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
        bool ignore_outer_gaps;
        
    private:
        int argc;
        char** argv;
};

#endif