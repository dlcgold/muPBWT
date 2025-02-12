//
// Created by dlcgold on 16/08/22.
//

#ifndef RLPBWT_COLUMN_NAIVE_MS_H
#define RLPBWT_COLUMN_NAIVE_MS_H

#include "utils.h"
#include <ostream>
#include <sdsl/int_vector.hpp>
#include <utility>
#include <vector>

class rl_column {
public:
  /**
   * @brief bool to check first value of the column in PBWT matrix
   * (assuming biallelic)
   */
  bool zero_first{};
  /**
   * @brief total number of zeros in the column
   */
  unsigned int count_0{};

  /**
   * @brief compressed sdsl int vector for positions of runs head
   */
  sdsl::int_vector<> p;

  /**
   * @brief compressed sdsl int vector for u/v value of a run
   */
  sdsl::int_vector<> uv;

  /**
   * @brief compressed sdsl int vector for thresholds
   */
  sdsl::int_vector<> t;

  /**
   * @brief sdsl compressed int vector for runs-head prefix array values
   */
  sdsl::int_vector<> sample_beg;

  /**
   * @brief sdsl compressed int vector for runs-tail prefix array values
   */
  sdsl::int_vector<> sample_end;
  /**
   * @brief sdsl compressed int vector for runs-head lcp array values
   */
  sdsl::int_vector<> sample_beg_lcp;

  sdsl::int_vector<> i_k;
  sdsl::int_vector<> l_k;
  sdsl::int_vector<> i_e_k;
  sdsl::int_vector<> l_e_k;
  /**
   * @brief constructor for a naive column supporting matching statistics
   * use
   *
   * @param zeroFirst
   * @param count0
   * @param p
   * @param uv
   * @param t
   * @param sample_beg
   * @param sample_end
   */
  rl_column(bool zeroFirst, unsigned int count0, sdsl::int_vector<> p,
            sdsl::int_vector<> uv, sdsl::int_vector<> t,
            sdsl::int_vector<> sample_beg, sdsl::int_vector<> sample_end,
            sdsl::int_vector<> sample_beg_lcp
            // ,sdsl::int_vector<> i_k,
            // sdsl::int_vector<> l_k, sdsl::int_vector<> i_e_k,
            // sdsl::int_vector<> l_e_k
  );

  /**
   * @brief default constructor
   */
  rl_column();

  /**
   * @brief default destructor
   */
  virtual ~rl_column();

  friend std::ostream &operator<<(std::ostream &os, const rl_column &naive);

  /**
   * @brief function to obtain size in bytes of the column
   * @param verbose bool for extra prints
   * @return size in bytes
   */
  unsigned long long size_in_bytes(bool verbose = false) const;

  /**
   * @brief function to obtain size in megabytes of the column
   * @param verbose bool for extra prints
   * @return size in megabytes
   */
  double size_in_mega_bytes(bool verbose = false) const;

  /**
   * @brief function to serialize the column object
   * @param out std::ostream object to stream the serialization
   * @return size of the serialization
   */
  size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                   const std::string &name = "");

  /**
   * @brief function to load the column object
   * @param in std::istream object from which load the column
   */
  void load(std::istream &in);
};

#endif // RLPBWT_COLUMN_NAIVE_MS_H
