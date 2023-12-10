#ifndef MS_K_H
#define MS_K_H

#include <iostream>
#include <ostream>
#include <sdsl/int_vector.hpp>
#include <tuple>
#include <utility>
#include <vector>

/**
 * @brief class to represent matching statistics
 */
class ms_k {
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
   * @brief support row vector of the matching statistics
   */
  std::vector<unsigned int> row_supp;

  /**
   * @brief support length vector of the matching statistics
   */
  std::vector<unsigned int> len_supp;

  /**
   * @brief default constructor
   */
  ms_k();

  /**
   * @brief default destructor
   */
  virtual ~ms_k();
  /**
   * @brief constructor from a couple of vectors
   * @param pos row vector of the matching statistics
   * @param len length vector of the matching statistics
   */
  ms_k(std::vector<unsigned int> row, std::vector<unsigned int> len,
       std::vector<unsigned int> row_supp, std::vector<unsigned int> len_supp);

  /**
   * @brief constructor from a number of matches
   * @param size number of matches
   */
  explicit ms_k(unsigned int size);

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

  friend std::ostream &operator<<(std::ostream &os, const ms_k &ms_k);
};

#endif // MS_K_H
