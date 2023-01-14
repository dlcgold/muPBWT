//
// Created by dlcgold on 20/02/22.
//

#ifndef RLPBWT_MS_MATCHES_H
#define RLPBWT_MS_MATCHES_H

#include <utility>
#include <iostream>
#include <vector>
#include <tuple>
#include <ostream>

/**
 * @brief class to represent matches obtained using matching statistics
 */
class ms_matches {
public:
    /**
     * @brief vector of tuple that represent match. In order:
     * - row obtained from row vector of the matching statistics
     * - length of the match obtained from length vector of the matching
     *   statistics
     * - ending column of the match
     */
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> basic_matches;

    /**
     * @brief vector that, for every match, represent the haplotypes that are
     * matching at that point
     */
    std::vector<std::vector<unsigned int>> haplos;

    /**
     * @brief default constructor
     */
    ms_matches();

    /**
     * @brief function to obtain size in bytes of the matches
     * @param verbose bool for extra prints
     * @return size in bytes
     */
    unsigned long long size_in_bytes() const;

    /**
     * @brief function to obtain size in megabytes of the matches
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes() const;
    friend std::ostream &operator<<(std::ostream &os, const ms_matches &ms);

};


#endif //RLPBWT_MS_MATCHES_H
