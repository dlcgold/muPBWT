//
// Created by dlcgold on 21/02/22.
//

#include "../include/matches_naive.h"


std::ostream &operator<<(std::ostream &os, const matches_naive &matches) {
    os << "\nmatches (<ending>, <length>, <number of haplotypes>):\n";
    for (const auto &basic_match: matches.basic_matches) {
        os << std::get<2>(basic_match) << ", " << std::get<1>(basic_match)
           << ", " << std::get<0>(basic_match);
        os << "\n";
    }
    return os;
}


unsigned long long matches_naive::size_in_bytes() const {
    unsigned long long size = 0;
    size += (sizeof(unsigned int) * this->basic_matches.size() * 3);
    return size;
}

double matches_naive::size_in_mega_bytes() const {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double) (sizeof(unsigned int) * this->basic_matches.size() * 3);
    return size * to_mega;
}

matches_naive::matches_naive() = default;
