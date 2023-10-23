#include <iostream>
#include <string>
#include "common.hpp"

void readargs(int argc, char *argv[], std::string &ini, std::string &circuit){
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg.length() > 1 && arg[0] == '-') {
            std::string option = arg.substr(1);

            // Handle different options
            if (option == "i") { // config file
                i++;
                std::cout << "[Config File]: " << argv[i] << std::endl;
                ini = std::string(argv[i]);
            } else if (option == "c") { // circuit file
                i++;
                std::cout << "[Circuit File]: " << argv[i] << std::endl;
                circuit = std::string(argv[i]);
            } else {
                // Handle unknown option
                std::cout << "[Unknown Option]: " << option << std::endl;
                exit(1);
            }
        } else {
            // Handle non-option argument
            std::cout << "[Non-option Argument]: " << arg << std::endl;
            exit(1);
        }
    }
}
