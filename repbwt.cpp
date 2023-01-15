#include <getopt.h>
#include <iostream>
#include "include/rlpbwt_int.h"

void printHelp() {
    std::cout << "Usage: PLOUTON [options]\n" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input_file <path>\t macs file for panel" << std::endl;
    std::cout << "  -s, --save <path>\t  path to save serialization "
              << std::endl;
    std::cout << "  -l, --load <path>\t path to load serialization"
              << std::endl;
    std::cout << "  -o, --output <path>\t path to query output" << std::endl;
    std::cout << "  -q, --query <path>\t path to macs query file" << std::endl;
//    std::cout << "  -e, --extend\t extend matches"
//              << std::endl;
    std::cout << "  -m, --macs\t use macs format for panel and queries"
              << std::endl;
    std::cout << "  -v, --verbose\t extra prints" << std::endl;
    std::cout << "  -V, --fverbose\t extra prints for functions (cautions)"
              << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

int main(int argc, char **argv) {
    // TODO option to print stats
    if (argc == 1) {
        printHelp();
        exit(EXIT_SUCCESS);
    }
    bool verbose = false;
    bool print_verbose = false;
    //bool extend = true;
    bool query = false;
    bool macs = false;
    std::string matrix_input = "";
    std::string memorize_file = "";
    std::string load_file = "";
    std::string output = "";
    std::string query_input = "";
    int c;

    while (true) {
        static struct option long_options[] = {
                {"input",    required_argument, nullptr, 'i'},
                {"save",     required_argument, nullptr, 's'},
                {"load",     required_argument, nullptr, 'l'},
                {"output",   required_argument, nullptr, 'o'},
                {"query",    required_argument, nullptr, 'q'},
                //{"extend",   no_argument,       nullptr, 'e'},
                {"macs",     no_argument,       nullptr, 'm'},
                {"fverbose", no_argument,       nullptr, 'V'},
                {"verbose",  no_argument,       nullptr, 'v'},
                {"help",     no_argument,       nullptr, 'h'},
                {nullptr,    0,                 nullptr, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "i:s:l:o:q:emvVh", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':
                matrix_input = optarg;
                break;
            case 's':
                memorize_file = optarg;
                break;
            case 'l':
                load_file = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'q':
                query_input = optarg;
                break;
            case 'm':
                macs = true;
                break;
            case 'v':
                print_verbose = true;
                break;
            case 'V':
                verbose = true;
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
            default:
                printHelp();
                exit(EXIT_FAILURE);
        }
    }

    if (!query_input.empty()) {
        query = true;
    }
    if (memorize_file.empty() && (query_input.empty() || output.empty())) {
        std::cerr << "Error: nothing to do\n";
        printHelp();
        exit(EXIT_FAILURE);
    }
    if (query && output.empty()) {
        std::cerr << "Error: output file required if query requested\n";
        printHelp();
        exit(EXIT_FAILURE);
    }


    if (matrix_input.empty() && load_file.empty()) {
        std::cerr << "Error: input or load file required\n";
        exit(EXIT_FAILURE);
    }

    if (query) {
        if (output.empty()) {
            std::cerr << "Error: output file required\n";
            exit(EXIT_FAILURE);
        }
        if (query_input.empty()) {
            std::cerr << "Error: query file required\n";
            exit(EXIT_FAILURE);
        }
    }


    clock_t START = clock();
    if (load_file.empty()) {

        rlpbwt_int rlpbwt(matrix_input.c_str(), verbose, macs);
        if (!memorize_file.empty()) {
            std::ofstream outstream;
            outstream.open(memorize_file.c_str());
            rlpbwt.serialize(outstream);
            outstream.close();
        }
        auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
        if (print_verbose) {
            std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                      << " bytes\n";
            std::cout
                    << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                    << " megabytes\n----\n";
            std::cout
                    << "estimated dense size: "
                    << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                    << " megabytes\n----\n";
        }
        std::cout << "built/loaded in: " << time_build << " s\n";
        if (query) {
            std::cout << "start querying \n";
            START = clock();
            rlpbwt.query_match(query_input.c_str(), output.c_str(), verbose,
                               macs);

            auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
            std::cout << "queried in " << time_query << " s\n";

        }
    } else {
        rlpbwt_int rlpbwt;
        std::ifstream load;
        load.open(load_file.c_str());
        rlpbwt.load(load);
        load.close();
        auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
        if (print_verbose) {
            std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                      << " bytes\n";
            std::cout
                    << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                    << " megabytes\n----\n";
            std::cout
                    << "estimated dense size: "
                    << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                    << " megabytes\n----\n";
        }
        std::cout << "built/loaded in: " << time_build << " s\n";
        if (query) {
            std::cout << "start querying \n";
            START = clock();
            rlpbwt.query_match(query_input.c_str(), output.c_str(), verbose,
                               macs);

            auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
            std::cout << "queried in " << time_query << " s\n";

        }
    }

    return 0;
}
