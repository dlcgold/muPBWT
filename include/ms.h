//
// Created by dlcgold on 19/02/22.
//

#ifndef RLPBWT_MS_H
#define RLPBWT_MS_H


#include <utility>
#include <iostream>
#include <vector>
#include <tuple>
#include <ostream>
#include <sdsl/int_vector.hpp>

/**
 * @brief class to represent matching statistics
 */
class ms {
public:
    /**
     * @brief row vector of the matching statistics
     */
    std::vector<unsigned int> row;

    /**
     * @brief length vector of the matching statistics
     */
    std::vector<unsigned int> len;

    /**
     * @brief default constructor
     */
    ms();

  /**
   * @brief default destructor
   */
  virtual ~ms();
    /**
     * @brief constructor from a couple of vectors
     * @param pos row vector of the matching statistics
     * @param len length vector of the matching statistics
     */
    ms(std::vector<unsigned int> pos, std::vector<unsigned int> len);

    /**
     * @brief constructor from a number of matches
     * @param size number of matches
     */
    explicit ms(unsigned int size);

    /**
     * @brief function to obtain size in bytes of the matching statistics
     * @param verbose bool for extra prints
     * @return size in bytes
     */
    unsigned long long size_in_bytes() const;

    /**
     * @brief function to obtain size in megabytes of the matching statistics
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes() const;

    friend std::ostream &operator<<(std::ostream &os, const ms &ms);
};


#endif //RLPBWT_MS_H
