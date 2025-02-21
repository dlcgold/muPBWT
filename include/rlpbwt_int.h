//
// Created by dlcgold on 16/08/22.
//

#ifndef RLPBWT_RLPBWT_NAIVE_MS_H
#define RLPBWT_RLPBWT_NAIVE_MS_H

#include "exceptions.h"
#include "htslib/vcf.h"
#include "ms.h"
#include "ms_k.h"
#include "ms_matches.h"
#include "phi_ds.h"
#include "rl_column.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <memory>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

std::vector<sdsl::sd_vector<>> tb(const std::vector<sdsl::bit_vector> &matrix) {
  if (matrix.empty())
    return {};

  size_t rows = matrix.size();
  size_t cols = matrix[0].size();

  std::vector<sdsl::bit_vector> transposed(cols, sdsl::bit_vector(rows, 0));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      transposed[j][i] = matrix[i][j];
    }
  }

  std::vector<sdsl::sd_vector<>> compressed_transposed;
  for (const auto &bv : transposed) {
    sdsl::sd_vector<> sdv(bv);
    compressed_transposed.push_back(sdv);
  }

  return compressed_transposed;
}

/**
 * @brief data structure for matching-statistics supported RLPBWT naive
 */
class rlpbwt_int {
public:
  /**
   * @brief bool to check if the RLPBWT is extended with phi/phi_inv structure
   */
  bool is_extended{};

  /**
   * @brief vector of matcing statistics supported naive columns
   */
  std::vector<rl_column> cols;

  /**
   * @brief phi/phi_inv support data structure
   */
  phi_ds *phi = nullptr;

  unsigned int k_smem = 0;

  std::vector<sdsl::sd_vector<>> panel;

private:
  /**
   * @brief compressed int vector for last prefix array
   */
  sdsl::int_vector<> last_pref;
  sdsl::int_vector<> last_div;

  rl_column build_column(std::string &column, std::vector<unsigned int> &pref,
                         std::vector<unsigned int> &div,
                         std::vector<intv> &supp_b, std::vector<intv> &supp_e,
                         unsigned int col_i = 0) {
    // unsigned int height = pref.size();
    //  variable for "c" value
    unsigned int count0 = 0;
    // variables for compute "u" and "v" values
    unsigned int u = 0;
    unsigned int v = 0;
    // temporary variables for compute "u" and "v" values
    unsigned int count0tmp = 0;
    unsigned int count1 = 0;
    // bool to check first symbol fo the column
    bool start = true;

    // update start and "c" value
    for (unsigned int i = 0; i < height; i++) {
      if (i == 0 && column[pref[i]] == '1') {
        start = false;
      }
      if (column[i] == '0') {
        count0++;
      }
    }

    // temporary variable to compute thresholds
    unsigned int lcs = 0;

    // initialize a vector of pair in order to build final sdsl int_vector
    // for p and u/v
    std::vector<std::pair<unsigned int, unsigned int>> rows;
    std::vector<unsigned int> thr;
    // support vector to store prefix array samples
    std::vector<std::pair<unsigned int, unsigned int>> samples;
    std::vector<unsigned int> samples_lcp;

    // temporary variable for p
    unsigned int p_tmp = 0;

    // temporary variable for prefix array value at the begin of a run
    unsigned int tmp_beg = 0;
    // temporary variable for divergence array value at the begin of a run
    unsigned int tmp_lcp = 0;
    // temporary variable for thresholds
    unsigned int tmp_thr = 0;

    // bools to check if we are at the beginning of a run and if we have to
    // swap the counting of zeros and ones
    bool begrun = true;
    bool pushz = false;
    bool pusho = false;

    // first swap according to first value of the column
    if (start) {
      pusho = true;
    } else {
      pushz = true;
    }

    // iteration over the entire column
    for (unsigned int i = 0; i < height; i++) {
      // if we are at the beginning of a run we save previous temporary
      // values for "u" and "v"
      if (begrun) {
        tmp_beg = pref[i];
        tmp_lcp = div[i];
        u = count0tmp;
        v = count1;
        begrun = false;
      }
      // increment temporary variables
      if (column[pref[i]] == '1') {
        count1++;
      } else {
        count0tmp++;
      }

      // updating thresholds (iff thresholds are required)

      if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
        tmp_thr = i;
        lcs = div[i];
      }
      if (div[i] < lcs) {
        tmp_thr = i;
        lcs = div[i];
      }

      // record starting position of a run
      if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
        p_tmp = i;
      }

      // do stuff at run change
      if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
        // update vector for p and u/v and swap the case to study in
        // next run
        if (pusho) {
          rows.emplace_back(p_tmp, v);
          std::swap(pusho, pushz);
        } else {
          rows.emplace_back(p_tmp, u);
          std::swap(pusho, pushz);
        }

        if (i + 1 != height && div[i + 1] < div[tmp_thr]) {
          thr.push_back(i + 1);
        } else {
          thr.push_back(tmp_thr);
        }

        samples.emplace_back(tmp_beg, pref[i]);
        samples_lcp.emplace_back(tmp_lcp);
        begrun = true;
      }
    }
    // create and compress sdsl int_vector

    sdsl::int_vector<> t_vec(rows.size());
    sdsl::int_vector<> p_vec(rows.size());
    sdsl::int_vector<> uv_vec(rows.size());
    sdsl::int_vector<> sb_vec(rows.size());
    sdsl::int_vector<> se_vec(rows.size());
    sdsl::int_vector<> sbl_vec(rows.size());
    sdsl::int_vector<> sel_vec(rows.size());
    for (unsigned int i = 0; i < rows.size(); i++) {
      p_vec[i] = rows[i].first;
      uv_vec[i] = rows[i].second;
      t_vec[i] = thr[i];
      sb_vec[i] = samples[i].first;
      se_vec[i] = samples[i].second;
      sbl_vec[i] = samples_lcp[i];
      supp_b[samples[i].first].push_back(col_i);
      supp_e[samples[i].second].push_back(col_i);
      // supp_l[samples_lcp[i]].push_back(count);
    }
    sdsl::util::bit_compress(p_vec);
    sdsl::util::bit_compress(uv_vec);
    sdsl::util::bit_compress(t_vec);
    sdsl::util::bit_compress(sb_vec);
    sdsl::util::bit_compress(se_vec);
    sdsl::util::bit_compress(sbl_vec);
    // return the column
    return {start, count0, p_vec, uv_vec, t_vec, sb_vec, se_vec, sbl_vec};
  }

  /**
   * @brief function to update prefix/divergence (lcp) arrays as in Durbin
   * @param column column of the panel
   * @param pref current prefix array
   * @param div current divergence array (as lcp array)
   */
  //    static void update(std::string &column, std::vector<unsigned int> &pref,
  //                       sdsl::int_vector<> &div) {
  static void update(std::string &column, std::vector<unsigned int> &pref,
                     std::vector<unsigned int> &div) {
    unsigned int height = pref.size();
    std::vector<unsigned int> new_pref(height);
    // sdsl::int_vector<> new_div(height);
    std::vector<unsigned int> new_div(height);
    unsigned int count0 = 0;
    unsigned int lcs = -1;

    for (unsigned int i = 0; i < height; i++) {
      // lcs = std::min(lcs, static_cast<unsigned int>(div[i]));
      lcs = std::min(lcs, div[i]);
      if (column[pref[i]] == '0') {
        new_pref[count0] = pref[i];
        new_div[count0] = lcs + 1;
        count0++;
        lcs = INT_MAX;
      }
    }

    int count1 = 0;
    lcs = -1;
    for (unsigned int i = 0; i < height; i++) {
      // lcs = std::min(lcs, static_cast<unsigned int>(div[i]));
      lcs = std::min(lcs, div[i]);
      if (column[pref[i]] == '1') {
        new_pref[count0 + count1] = pref[i];
        new_div[count0 + count1] = lcs + 1;
        count1++;
        lcs = INT_MAX;
      }
    }
    new_div[0] = 0;
    if (count0 != height) {
      new_div[count0] = 0;
    }
    div = new_div;
    pref = new_pref;
  }

  /**
   * @brief function to find the subvector os size k which contain the maximum
   * minimum value
   * @param vector input vector
   * @param k k-value
   */
  std::pair<unsigned int, int>
  maxSubvector(const std::vector<unsigned int> &vector, int k) {

    int size = vector.size();
    // if (k > vector.size()) {
    //   k = vector.size();
    // }

    int val = std::numeric_limits<int>::min();
    int ind = -1;

    std::deque<int> deq;

    // using a sliding window to keep track of all the subvectors
    for (int i = 0; i < size; ++i) {
      // Remove elements outside the current window
      while (!deq.empty() && deq.front() < i - k + 1) {
        deq.pop_front();
      }

      // Remove elements in the current window that are greater than the current
      // element
      while (!deq.empty() && vector[deq.back()] > vector[i]) {
        deq.pop_back();
      }

      // Add the current element index to the window
      deq.push_back(i);

      // Check if the window is of size k
      if (i >= k - 1) {
        // The front of the deque contains the index of the minimum element in
        // the current window
        int min_value = vector[deq.front()];

        // Update if the current minimum is larger than the previous maximum
        // minimum
        if (min_value > val) {
          val = min_value;
          ind = i;
        }
      }
    }

    return {val, ind};
  }

  /**
   * @brief function to computer k intervals anbd its min length at run boundary
   * @param div divergence array
   * @param verbose bool for verbose print
   */
  void update_k_intervals(const std::vector<unsigned int> &div,
                          bool verbose = false) {
    auto h = div.size();

    // iterate over every run boundary
    for (unsigned int i = 0; i < this->cols.back().p.size(); i++) {
      // run-head

      // perform a lf-mapping to the next column
      auto lf_pos =
          rlpbwt_int::lf(this->cols.size() - 1, this->cols.back().p[i],
                         get_next_char(this->cols.back().zero_first,
                                       index_to_run(this->cols.back().p[i],
                                                    this->cols.size() - 1)));

      // get the interval [lf(r)-k-2,lf(r)+k-1] eventually fixing at panel
      // boundaries
      auto a = 0;
      if (lf_pos + 2 >= this->k_smem) {
        a = lf_pos - this->k_smem + 2;
      }
      auto e = h - 1;
      if (lf_pos - 1 + this->k_smem < h) {
        e = lf_pos + this->k_smem - 1;
      }
      // get the subvector of div and compute the subsubvector os size k-1 with
      // max min value
      std::vector<unsigned int> subvector(div.begin() + a, div.begin() + e + 1);
      auto sub_m = maxSubvector(subvector, this->k_smem - 1);
      if (verbose) {
        std::cout << "a " << a << " e " << e << " sub_m.first " << sub_m.first
                  << " sub_m.second " << sub_m.second << "\n";
      }
      // sentiler is subvector short than k-1
      if (e - a + 1 < this->k_smem - 1) {
        sub_m.second = -1;
      }

      // compute the offset between the starting point of the subvector and lf
      // and the len
      if (sub_m.second != -1) {
        sub_m.second = (sub_m.second + a) - (k_smem - 2);
        if (sub_m.second > 0) {
          sub_m.second--;
        }
        this->cols.back().i_k[i] =
            static_cast<unsigned int>(sub_m.second) >= lf_pos
                ? 0
                : lf_pos - sub_m.second;

        this->cols.back().l_k[i] = sub_m.first;
      } else {
        this->cols.back().i_k[i] = 0;
        this->cols.back().l_k[i] = 0;
      }

      if (verbose) {
        std::cout << "col s " << this->cols.size() - 1 << " p "
                  << this->cols.back().p[i] << " lf_pos " << lf_pos << " a "
                  << a << " e " << e << " sub_m " << sub_m.first
                  << " sub_m.second " << sub_m.second << " i_k "
                  << this->cols.back().i_k[i] << " l_k "
                  << this->cols.back().l_k[i] << " using "
                  << this->cols.back().p[i] << " for "
                  << get_next_char(this->cols.back().zero_first,
                                   index_to_run(this->cols.back().p[i],
                                                this->cols.size() - 1))
                  << " -> ";
        for (unsigned int j = 0; j < subvector.size(); j++) {
          std::cout << subvector[j] << " ";
        }
        std::cout << "\n";
      }

      // run-tail

      // get lf checing if we are atthe end of the panel the do as before
      lf_pos = 0;
      if (i != this->cols.back().p.size() - 1) {
        lf_pos = rlpbwt_int::lf(
            this->cols.size() - 1, this->cols.back().p[i + 1] - 1,
            get_next_char(this->cols.back().zero_first,
                          index_to_run(this->cols.back().p[i + 1] - 1,
                                       this->cols.size() - 1)));
      } else {
        lf_pos =
            rlpbwt_int::lf(this->cols.size() - 1, this->height - 1,
                           get_next_char(this->cols.back().zero_first,
                                         index_to_run(this->height - 1,
                                                      this->cols.size() - 1)));
      }

      a = 0;
      if (lf_pos + 2 >= this->k_smem) {
        a = lf_pos - this->k_smem + 2;
      }
      e = h - 1;
      if (lf_pos - 1 + this->k_smem < h) {
        e = lf_pos + this->k_smem - 1;
      }

      std::vector<unsigned int> subvector2(div.begin() + a,
                                           div.begin() + e + 1);
      sub_m = maxSubvector(subvector2, this->k_smem - 1);
      if (verbose) {
        std::cout << "a " << a << " e " << e << " sub_m.first " << sub_m.first
                  << " sub_m.second " << sub_m.second << "\n";
      }

      if (e - a + 1 < this->k_smem - 1) {
        sub_m.second = -1;
      }
      if (sub_m.second != -1) {
        sub_m.second = (sub_m.second + a) - (k_smem - 2);
        if (sub_m.second > 0) {
          sub_m.second--;
        }
        this->cols.back().i_e_k[i] =
            static_cast<unsigned int>(sub_m.second) >= lf_pos
                ? 0
                : lf_pos - sub_m.second;

        this->cols.back().l_e_k[i] = sub_m.first;
      } else {
        this->cols.back().i_e_k[i] = 0;
        this->cols.back().l_e_k[i] = 0;
      }
      if (verbose) {
        std::cout << "col e " << this->cols.size() - 1 << " p "
                  << this->cols.back().p[i] << " lf_pos " << lf_pos << " a "
                  << a << " e " << e << " sub_m " << sub_m.first
                  << " sub_m.second " << sub_m.second << " i_e_k "
                  << this->cols.back().i_e_k[i] << " l_e_k "
                  << this->cols.back().l_e_k[i] << " -> ";
        for (unsigned int j = 0; j < subvector2.size(); j++) {
          std::cout << subvector2[j] << " ";
        }
        std::cout << "\n";
      }
    }

    // bit-compress all the values
    sdsl::util::bit_compress(this->cols.back().i_k);
    sdsl::util::bit_compress(this->cols.back().l_k);
    sdsl::util::bit_compress(this->cols.back().i_e_k);
    sdsl::util::bit_compress(this->cols.back().l_e_k);
  }
  /**
   * @brief function to compute the lf mapping, w(i, s) function in Durbin
   * @param col_index index of the column
   * @param index virtual index of the row of the original panel
   * @param symbol symbol s
   * @param verbose bool for extra print
   * @return the index computed with the lf-mapping
   */
  unsigned int lf(unsigned int col_index, unsigned int index, char symbol,
                  bool verbose = false) const {
    // get run
    unsigned int run_index = index_to_run(index, col_index);
    // get offset
    unsigned int offset = index - this->cols[col_index].p[run_index];
    // undoing the offsets when they are wrong/useless
    if ((symbol == '0' &&
         get_next_char(this->cols[col_index].zero_first, run_index) == '1') ||
        (symbol == '1' &&
         get_next_char(this->cols[col_index].zero_first, run_index) == '0')) {
      offset = 0;
    }
    // obtain "u" and "v"
    auto uv = uvtrick(col_index, run_index);
    if (verbose) {
      std::cout << uv.first << ", " << uv.second << "\n";
    }
    // fix for the last index that's "outside" the column
    if (this->cols[col_index].p[run_index] + offset == this->height) {
      if (get_next_char(this->cols[col_index].zero_first, run_index) == '0') {
        uv.second--;
      } else {
        uv.first--;
      }
    }
    // computer lf-mapping as Durbin's w(i, s), eventually using offset
    if (symbol == '0') {
      return uv.first + offset;
    } else {
      return this->cols[col_index].count_0 + uv.second + offset;
    }
  }

  unsigned int reverse_lf(unsigned int col_index, unsigned int index,
                          bool verbose) const {
    // TODO implement as binary search
    // by design if we try to work on first column the function return 0
    if (col_index == 0) {
      return 0;
    }
    // we extract the "c" value from the previous column
    col_index = col_index - 1;
    unsigned int c = this->cols[col_index].count_0;
    // initialize u/v, offset and the new run index
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int offset = 0;
    unsigned int pos = 0;

    // bool to interrupt the search
    bool found = false;
    if (verbose) {
      std::cout << "c: " << c << "\n";
    }
    // two cases:
    // - if index is less than previous "c" it means that it comes from a
    // zero
    //   element, so we will search the correct u value
    // - otherwise it means that it comes from a one element, so we will
    // search
    //   the correct v value
    if (index < c) {
      u = index;
      if (verbose) {
        std::cout << "u: " << u << "\n";
      }
      // iteration over the u values to find the correct one
      int s = 0;
      int e = this->cols[col_index].p.size() - 2;
      int m = 0;
      unsigned int prevu = 0;
      unsigned int nextu = 0;
      while (s <= e) {
        m = (s + e) / 2;
        prevu = uvtrick(col_index, m).first;
        nextu = uvtrick(col_index, m + 1).first;
        if (prevu <= u && u < nextu) {
          pos = m;
          found = true;
          break;
        } else if (prevu <= u) {
          s = m + 1;
        } else {
          e = m - 1;
        }
      }

      // if not found we are at the last run
      if (!found) {
        pos = this->cols[col_index].p.size() - 1;
      }
      if (verbose) {
        std::cout << "row: " << pos << "\n";
      }

      // using the correct u value in previous column to obtain the previous
      // index using offset
      unsigned int curru = uvtrick(col_index, pos).first;
      offset = u - curru;
      if (verbose) {
        std::cout << "offset: " << offset << "\n";
      }
      return this->cols[col_index].p[pos] + offset;
    } else {
      int s = 0;
      int e = this->cols[col_index].p.size() - 2;
      int m = 0;
      v = index - c;
      if (verbose) {
        std::cout << "v: " << v << "\n";
      }
      unsigned int prevv = 0;
      unsigned int nextv = 0;
      while (s <= e) {
        m = (s + e) / 2;
        prevv = uvtrick(col_index, m).second;
        nextv = uvtrick(col_index, m + 1).second;
        if (prevv <= v && v < nextv) {
          pos = m;
          found = true;
          break;
        } else if (prevv <= v) {
          s = m + 1;
        } else {
          e = m - 1;
        }
      }
      if (verbose) {
        std::cout << "s " << s << " e " << e << " m " << m << " pos " << pos
                  << "\n";
      }

      // if not found we are at the last run
      if (!found) {
        pos = this->cols[col_index].p.size() - 1;
      }
      if (verbose) {
        std::cout << "row: " << pos << "\n";
      }
      // using the correct v value in previous column to obtain the previous
      // index using offset
      unsigned int currv = uvtrick(col_index, pos).second;
      offset = v - currv;
      if (verbose) {
        std::cout << "offset: " << offset << "\n";
      }
      return this->cols[col_index].p[pos] + offset;
    }
    return 0;
  }
  /**
   * @brief function to map an index to the correct run in a column
   * @param index index to map
   * @param col_index column index
   * @return run index
   */
  unsigned int index_to_run(unsigned int index, unsigned int col_index) const {
    // if requested index is equal or greater than the p value of the last
    // run return the index of the last run
    if (index >= this->cols[col_index].p[this->cols[col_index].p.size() - 1]) {
      return this->cols[col_index].p.size() - 1;
    }

    // binary search to compute run index
    unsigned int bi = 0;
    unsigned int e = this->cols[col_index].p.size();
    unsigned int pos = (e - bi) / 2;
    while (pos != e && this->cols[col_index].p[pos] != index) {
      if (index < (unsigned int)this->cols[col_index].p[pos]) {
        e = pos;
      } else {
        if (pos + 1 == e ||
            (unsigned int)this->cols[col_index].p[pos + 1] > index) {
          break;
        }
        bi = pos + 1;
      }
      pos = bi + (e - bi) / 2;
    }
    return pos;
  }

  /**
   * @brief trick to extract u and v value from a run in rlpbwt column
   * @param col_index index of the column
   * @param index virtual index of the row of the original panel
   * @return a std::pair with u as first and v as second
   */
  std::pair<unsigned int, unsigned int> uvtrick(unsigned int col_index,
                                                unsigned int run_index) const {
    unsigned int u;
    unsigned int v;
    // if run index is 0 u = v = 0
    // in other case, based on first symbol of the column
    // we have u/v in the same row and v/u in the previous one
    if (run_index == 0) {
      u = 0;
      v = 0;
    } else if (run_index % 2 == 0) {
      u = this->cols[col_index].uv[run_index - 1];
      v = this->cols[col_index].uv[run_index];
      if (!this->cols[col_index].zero_first) {
        std::swap(u, v);
      }
    } else {
      u = this->cols[col_index].uv[run_index];
      v = this->cols[col_index].uv[run_index - 1];
      if (!this->cols[col_index].zero_first) {
        std::swap(u, v);
      }
    }
    return {u, v};
  }

  /**
   * @brief function to extend matching statistics matches with matching
   * rows indices
   * @param ms_matches matching statistics matches that will be extended
   */
  void extend_haplos(ms_matches &ms_matches,
                     std::vector<unsigned int> ms_supp) {

    // iterate over every basic match
    for (unsigned int i = 0; i < ms_matches.basic_matches.size(); i++) {
      // initialize the vector that will contain row indices
      std::vector<unsigned int> haplos;

      // extract information from the current basic match
      unsigned int start_row = std::get<0>(ms_matches.basic_matches[i]);
      haplos.emplace_back(start_row);
      unsigned int curr_len = std::get<1>(ms_matches.basic_matches[i]);
      unsigned int curr_col = std::get<2>(ms_matches.basic_matches[i]);

      if (curr_len == 1) {

        char curr_s = get_next_char(this->cols[curr_col].zero_first,
                                    index_to_run(ms_supp[curr_col], curr_col));
        auto full_col = this->get_col(curr_col);
        auto pa = this->get_prefix(curr_col);

        for (unsigned int c = 0; c < full_col.size(); c++) {
          if (full_col[c] == curr_s)
            haplos.emplace_back(pa[c]);
        }
        // std::cout << curr_col << "\n";
        ms_matches.haplos.emplace_back(haplos);
        continue;
      }
      // initialize boolean and temporary variables for go up/down in
      // search of matching rows
      bool check_down = true;
      unsigned int down_row = 0;
      bool check_up = true;
      unsigned int up_row = 0;

      char curr_s = get_next_char(this->cols[curr_col].zero_first,
                                  index_to_run(ms_supp[curr_col], curr_col));
      unsigned int down_index = lf(curr_col, ms_supp[curr_col], curr_s) + 1;
      unsigned int up_index = lf(curr_col, ms_supp[curr_col], curr_s);
      // go down/up and add row that has a lce of at least current len in
      // matching statistics (if panel is not saved as SLP the computation
      // of lce is simulated with random access to the panel)

      // down
      while (check_down) {
        auto phi_res = this->phi->phi_inv(start_row, curr_col + 1);
        if (!phi_res.has_value()) {
          break;
        }
        down_row = phi_res.value();
        if (this->phi->plcp(down_row, curr_col + 1) >= curr_len) {
          haplos.emplace_back(down_row);
          start_row = down_row;
        } else {
          check_down = false;
        }
        down_index++;
      }

      start_row = std::get<0>(ms_matches.basic_matches[i]);
      // up
      while (check_up) {
        auto phi_res = this->phi->phi(start_row, curr_col + 1);

        if (!phi_res.has_value()) {
          break;
        }
        up_row = phi_res.value();
        if (this->phi->plcp(start_row, curr_col + 1) >= curr_len) {
          haplos.emplace_back(up_row);
          start_row = up_row;
        } else {
          check_up = false;
        }
        up_index--;
      }
      // record the haplotypes
      ms_matches.haplos.emplace_back(haplos);
    }
  }

  std::vector<unsigned int> extend_single(unsigned int start_row,
                                          unsigned int curr_len,
                                          unsigned int curr_col,
                                          unsigned int curr_sup) {
    std::vector<unsigned int> haplos;
    if (curr_len == 0) {
      return haplos;
    }
    auto s_tmp = start_row;
    // extract information from the current basic match
    haplos.emplace_back(start_row);
    // initialize boolean and temporary variables for go up/down in
    // search of matching rows
    bool check_down = true;
    unsigned int down_row = 0;
    bool check_up = true;
    unsigned int up_row = 0;

    char curr_s = get_next_char(this->cols[curr_col].zero_first,
                                index_to_run(curr_sup, curr_col));
    unsigned int down_index = lf(curr_col, curr_sup, curr_s) + 1;
    unsigned int up_index = lf(curr_col, curr_sup, curr_s);
    // go down/up and add row that has a lce of at least current len in
    // matching statistics (if panel is not saved as SLP the computation
    // of lce is simulated with random access to the panel)

    // down
    while (check_down) {
      auto phi_res = this->phi->phi_inv(start_row, curr_col + 1);
      if (!phi_res.has_value()) {
        break;
      }
      down_row = phi_res.value();
      if (this->phi->plcp(down_row, curr_col + 1) >= curr_len) {
        haplos.emplace_back(down_row);
        start_row = down_row;
      } else {
        check_down = false;
      }
      down_index++;
    }

    start_row = s_tmp;
    // up
    while (check_up) {
      auto phi_res = this->phi->phi(start_row, curr_col + 1);

      if (!phi_res.has_value()) {
        break;
      }
      up_row = phi_res.value();
      if (this->phi->plcp(start_row, curr_col + 1) >= curr_len) {
        haplos.emplace_back(up_row);
        start_row = up_row;
      } else {
        check_up = false;
      }
      up_index--;
    }
    // record the haplotypes
    return haplos;
  }

  unsigned int lce(unsigned int q, unsigned int s, unsigned int c,
                   unsigned int l) {
    int col = static_cast<int>(c);
    unsigned int lce = 0;
    while (
        col >= 0 &&
        get_next_char(this->cols[col].zero_first, index_to_run(q, col)) ==
            get_next_char(this->cols[col].zero_first, index_to_run(s, col))) {
      if (col > 0) {
        q = reverse_lf(col, q, false);
        s = reverse_lf(col, s, false);
      }
      col--;
      lce++;
      if (lce == l) {
        return l;
      }
    }
    return lce;
  }

public:
  /**
   * @brief height of the panel
   */
  unsigned int width{};

  /**
   * @brief width of the panel
   */
  unsigned int height{};

  /**
   * @brief default constructor
   */
  rlpbwt_int() = default;

  /**
   * @brief default destructor
   */
  ~rlpbwt_int() {
    if (this->is_extended) {
      delete phi;
    }
  }

  /**
   * @brief constructor of a RLPBWT that support matching statistics
   * @param filename file with the panel
   * @param thr bool to enable thresholds computation
   * @param verbose bool fro extra prints
   */
  /**
   * @brief constructor of a RLPBWT that support matching statistics
   * @param filename file with the panel
   * @param thr bool to enable thresholds computation
   * @param verbose bool fro extra prints
   */
  explicit rlpbwt_int(const char *filename, std::unordered_set<unsigned int> ps,
                      bool verbose = false, int k_smem = 2) {
    this->k_smem = k_smem;
    std::vector<sdsl::bit_vector> panel_tmp;
    htsFile *fp = hts_open(filename, "rb");
    std::cout << "Reading VCF file...\n";
    if (fp == NULL) {
      throw FileNotFoundException{};
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();
    this->height = bcf_hdr_nsamples(hdr) * 2;
    this->width = 0;
    // while (bcf_read(fp, hdr, rec) >= 0) {
    //     this->width++;
    // }
    std::cout << "Total haplotypes: " << this->height << "\n";
    // std::cout << "h: " << this->height << "\n";
    // std::cout << "w: " << this->width << "\n";
    //  bcf_hdr_destroy(hdr);
    // hts_close(fp);
    //  bcf_destroy(rec);
    fp = hts_open(filename, "rb");
    hdr = bcf_hdr_read(fp);
    rec = bcf_init();

    // this->cols = std::vector<rl_column>(this->width + 1);
    std::vector<unsigned int> pref(this->height);
    std::vector<unsigned int> div(this->height);
    auto supp_b = std::vector<intv>(this->height);
    auto supp_e = std::vector<intv>(this->height);
    //            sdsl::int_vector<>
    //                    div(this->height);

    this->last_pref.resize(this->height);
    this->last_div.resize(this->height);

    for (unsigned int i = 0; i < this->height; i++) {
      pref[i] = i;
      div[i] = 0;
    }

    unsigned int count = 0;
    std::string last_col;
    std::string last_column;
    std::string new_column;

    // iterate each vcf record
    while (bcf_read(fp, hdr, rec) >= 0) {
      // if (auto search = ps.find(rec->pos); search == ps.end())
      //   continue;
      this->width++;
      // std::cout << rec->pos << "\n";
      std::cout << "Analyze site: " << count << "\r";
      new_column = "";
      sdsl::bit_vector tmp(this->height);
      int bi = 0;
      bcf_unpack(rec, BCF_UN_ALL);
      // read SAMPLE
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
      ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
      int max_ploidy = ngt / nsmpl;
      for (i = 0; i < nsmpl; i++) {
        int32_t *ptr = gt_arr + i * max_ploidy;
        for (j = 0; j < max_ploidy; j++) {
          // if true, the sample has smaller ploidy
          if (ptr[j] == bcf_int32_vector_end)
            break;

          // missing allele
          if (bcf_gt_is_missing(ptr[j]))
            exit(-1);

          // the VCF 0-based allele index
          int allele_index = bcf_gt_allele(ptr[j]);
          new_column += std::to_string(allele_index);
          tmp[bi] = allele_index;
          bi++;
        }
      }
      free(gt_arr);
      auto col = rlpbwt_int::build_column(new_column, pref, div, supp_b, supp_e,
                                          count);
      panel_tmp.push_back(tmp);
      // this->cols[count] = col;

      this->cols.emplace_back(col);
      rlpbwt_int::update(new_column, pref, div);
      if (k_smem <= 1) {
        this->cols.back().i_k.resize(0);
        this->cols.back().l_k.resize(0);
        this->cols.back().i_e_k.resize(0);
        this->cols.back().l_e_k.resize(0);
      } else if (k_smem > 1 && count >= 1) {
        rlpbwt_int::update_k_intervals(div);
      }
      last_col = new_column;
      count++;
    }
    this->panel = tb(panel_tmp);
    std::cout << std::endl;
    std::cout << "Total sites: " << this->width << "\n";
    for (unsigned int i = 0; i < pref.size(); i++) {
      this->last_pref[i] = pref[i];
      this->last_div[i] = div[i];
    }
    auto col =
        rlpbwt_int::build_column(last_col, pref, div, supp_b, supp_e, count);
    this->cols.emplace_back(col);
    for (unsigned int i = 0; i < this->height; i++) {
      if (supp_b[i].v.size() == 0 ||
          supp_b[i].v[supp_b[i].size - 1] != this->width) {
        supp_b[i].push_back(this->width);
      }
      supp_b[i].compress();
      if (supp_e[i].v.size() == 0 ||
          supp_e[i].v[supp_e[i].size - 1] != this->width) {
        supp_e[i].push_back(this->width);
      }
      supp_e[i].compress();
    }
    sdsl::util::bit_compress(this->last_pref);
    sdsl::util::bit_compress(this->last_div);
    this->phi =
        new phi_ds(this->cols, this->height, this->width, this->last_pref,
                   this->last_div, supp_b, supp_e, verbose);
    this->is_extended = true;

    bcf_hdr_destroy(hdr);
    hts_close(fp);
    bcf_destroy(rec);
  }

  /**
   * @brief function to delete the phi/phi_inv structure
   */
  void unextend() {
    if (this->is_extended) {
      this->phi = nullptr;
      this->is_extended = false;
    }
  }

  std::vector<unsigned int> slice_sd_vector(const sdsl::sd_vector<> &sdv,
                                            size_t X, size_t Y) {
    if (X > Y || Y >= sdv.size()) {
      std::cerr << "Invalid slice range!" << std::endl;
      std::cerr << X << " " << Y << " " << sdv.size() << std::endl;
      return {};
    }

    sdsl::sd_vector<>::rank_1_type rank(&sdv);
    sdsl::sd_vector<>::select_1_type select(&sdv);

    size_t rank_X = rank(X);
    // size_t rank_Y = rank(Y + 1);
    size_t rank_Y = (Y == sdv.size() - 1) ? rank(Y) + sdv[Y] : rank(Y + 1);
    std::vector<unsigned int> slice(Y - X + 1, 0);

    for (size_t i = rank_X; i < rank_Y; ++i) {
      size_t pos = select(i + 1);
      if (pos >= X && pos <= Y) {
        slice[pos - X] = 1;
      }
    }

    return slice;
  }
  /**
   * @brief function to compute matching statistics matches with a given
   * query using thresholds
   * @param query haplotype query as std::string
   * @param extend_matches bool to check if extend matching statistics
   * matches with rows
   * @param verbose bool for extra prints
   * @return matching statistics matches
   */
  ms_matches match_thr(const std::string &query, bool extend_matches = false,
                       bool verbose = false) {

    // compute the match iff |query| is equal to the width of the panel
    if (query.size() != this->width) {
      std::cout << query.size() << " != " << this->width << "\n";
      throw NotEqualLengthException{};
    }

    if (this->k_smem <= 1) {
      // for (auto c : this->cols) {
      //   std::cout << c;
      // }
      // initialize matching statistics
      auto ms_tot = compute_ms(query, extend_matches, verbose);
      auto ms = ms_tot.first;
      auto ms_supp = ms_tot.second;
      // initialize struct for matches
      ms_matches ms_matches;
      // save every match from matching statistics (when we have a "peak" in
      // ms len vector)
      for (unsigned int i = 0; i < ms.len.size(); i++) {
        if ((i != ms.len.size() - 1 && ms.len[i] > 0 &&
             ms.len[i] >= ms.len[i + 1]) ||
            (i == ms.len.size() - 1 && ms.len[i] != 0)) {
          ms_matches.basic_matches.emplace_back(ms.row[i], ms.len[i], i);
        }
      }
      // compute every row that are matching if required
      if (extend_matches) {
        if (verbose) {
          std::cout << "\nextending\n";
        }
        extend_haplos(ms_matches, ms_supp);
      }
      if (verbose) {
        std::cout << ms << "\n";
        std::cout << ms_matches << "\n";
      }

      return ms_matches;
    } else {
      // initialize matching statistics
      auto ms_tot = compute_ms_k(query, extend_matches, verbose);
      // for (auto c : this->cols) {
      //   std::cout << c;
      // }

      auto ms = ms_tot.first;
      auto ms_supp = ms_tot.second;
      if (verbose) {
        std::cout << ms_tot.first;
        std::cout << "\nms_supp:\t";
        for (auto m : ms_supp) {
          std::cout << m << "\t";
        }
        std::cout << "\n";
      }
      // initialize struct for matches
      ms_matches ms_matches;
      std::vector<unsigned int> f_len(ms.len.size());
      for (unsigned int i = 0; i < ms.len.size(); i++) {
        f_len[i] = std::min(ms.len[i], ms.len_supp[i]);
      }
      for (unsigned int i = 0; i < f_len.size(); i++) {
        if ((i != f_len.size() - 1 && f_len[i] > 0 &&
             f_len[i] >= f_len[i + 1]) ||
            (i == f_len.size() - 1 && f_len[i] != 0)) {
          ms_matches.basic_matches.emplace_back(ms.row[i], f_len[i], i);
        }
      }
      // save every match from matching statistics (when we have a "peak" in
      // ms len vector)
      // for (unsigned int i = 0; i < ms.len.size(); i++) {
      //   auto len_i = std::min(ms.len[i], ms.len_supp[i]);
      //   auto len_j = 0;
      //   if (i != ms.len.size() - 1) {
      //     len_j = std::min(ms.len[i + 1], ms.len_supp[i + 1]);
      //   }
      //   if ((i != ms.len.size() - 1 && len_i > 0 && len_i >= len_j) ||
      //       (i == ms.len.size() - 1 && len_i != 0)) {
      //     ms_matches.basic_matches.emplace_back(ms.row[i], len_i, i);
      //   }
      // }
      // compute every row that are matching if required
      if (extend_matches) {
        if (verbose) {
          std::cout << "\nextending\n";
        }
        extend_haplos(ms_matches, ms_supp);
      }
      auto i = 0;
      for (auto m : ms_matches.haplos) {
        if (m.size() < this->k_smem) {
          std::cerr << "ERROR: " << m.size() << " < " << this->k_smem << "\n";
          // std::cerr << ms_matches << "\n";
          for (auto e : m) {
            std::cerr << e << " ";
          }
          std::cerr << std::endl;
          std::cerr << std::get<0>(ms_matches.basic_matches[i]) << " "
                    << std::get<2>(ms_matches.basic_matches[i]) << " "
                    << std::get<1>(ms_matches.basic_matches[i]) << "\n";
          std::cerr << std::endl;
          // break;
          exit(1);
        }
        assert(m.size() >= this->k_smem);
        i++;
      }
      if (verbose) {
        // std::cout << ms << "\n";
        std::cout << ms_matches << "\n";
      }

      return ms_matches;
    }
  }

  ms_matches left_mpsc(const std::string &query, bool extend_matches = false,
                       bool verbose = false) {
    if (query.size() != this->width) {
      std::cout << query.size() << " != " << this->width << "\n";
      throw NotEqualLengthException{};
    }

    if (this->k_smem != 1) {
      // verbose = true;
      auto ms_tot = compute_ms_k(query, extend_matches, verbose);

      auto ms = ms_tot.first;
      // std::cout << ms << std::endl;
      auto ms_supp = ms_tot.second;
      // initialize struct for matches
      ms_matches ms_matches;
      // save every match from matching statistics (when we have a "peak" in
      // ms len vector)
      int j = this->width - 1;
      int jp = 0;
      std::vector<unsigned int> f_len(ms.len.size());
      for (unsigned int i = 0; i < ms.len.size(); i++) {
        f_len[i] = std::min(ms.len[i], ms.len_supp[i]);
      }
      while (j >= 0) {
        // std::cout << "at: "<< j << "  " << jp << " " << ms.len[j] <<
        // std::endl;
        if (f_len[j] != 0) {
          jp = j - f_len[j];
          // std::cout << "new jp: " << jp << "with j: " << j << std::endl;
          ms_matches.basic_matches.emplace_back(ms.row[j], j - (jp + 1) + 1, j);
          // std::cout << "add: "<< ms.row[j] << "  " << j - (jp + 1) + 1 << "
          //"
          //  << j << std::endl;
          j = jp;
        } else {
          jp = j - 1;
          j = jp;
        }
      }

      //  compute every row that are matching if required
      if (extend_matches) {
        if (verbose) {
          std::cout << "\nextending\n";
        }
        extend_haplos(ms_matches, ms_supp);
      }
      if (verbose) {
        std::cout << ms << "\n";
        std::cout << ms_matches << "\n";
      }

      return ms_matches;
    } else {
      // verbose = true;
      auto ms_tot = compute_ms(query, extend_matches, verbose);

      auto ms = ms_tot.first;
      std::cout << ms << std::endl;
      auto ms_supp = ms_tot.second;
      // initialize struct for matches
      ms_matches ms_matches;
      // save every match from matching statistics (when we have a "peak" in
      // ms len vector)
      int j = this->width - 1;
      int jp = 0;
      std::vector<unsigned int> f_len(ms.len.size());
      for (unsigned int i = 0; i < ms.len.size(); i++) {
        f_len[i] = ms.len[i];
      }
      while (j >= 0) {
        // std::cout << "at: "<< j << "  " << jp << " " << ms.len[j] <<
        // std::endl;
        if (f_len[j] != 0) {
          jp = j - f_len[j];
          // std::cout << "new jp: " << jp << "with j: " << j << std::endl;
          ms_matches.basic_matches.emplace_back(ms.row[j], j - (jp + 1) + 1, j);
          // std::cout << "add: "<< ms.row[j] << "  " << j - (jp + 1) + 1 << "
          //"
          //  << j << std::endl;
          j = jp;
        } else {
          jp = j - 1;
          j = jp;
        }
      }

      //  compute every row that are matching if required
      if (extend_matches) {
        if (verbose) {
          std::cout << "\nextending\n";
        }
        extend_haplos(ms_matches, ms_supp);
      }
      if (verbose) {
        std::cout << ms << "\n";
        std::cout << ms_matches << "\n";
      }

      return ms_matches;
    }
  }

  /**
   * @brief function to compute longer run in a column with a certain symbol
   * @param col_index column index
   * @param symbol run symbol
   */
  unsigned int best_run(unsigned int col_index, char symbol) {
    auto best_run_index = 0;
    if (this->cols[col_index].p.size() == 1) {
      return 0;
    }
    if (this->cols[col_index].p.size() == 2) {
      if ((symbol == '0' && this->cols[col_index].zero_first) ||
          (symbol == '1' && !this->cols[col_index].zero_first)) {
        return 0;
      } else {
        return 1;
      }
    }
    auto best_run_len = 0;
    if ((symbol == '0' && this->cols[col_index].zero_first) ||
        (symbol == '1' && !this->cols[col_index].zero_first)) {
      for (unsigned int i = 0; i < this->cols[col_index].p.size(); i += 2) {
        auto len_tmp = 0;
        if (i == this->cols[col_index].p.size() - 1) {
          len_tmp = this->height - this->cols[col_index].p[i];
        } else {
          len_tmp = this->cols[col_index].p[i + 1] - this->cols[col_index].p[i];
        }
        if (len_tmp > best_run_len) {
          best_run_len = len_tmp;
          best_run_index = i;
        }
      }
      return best_run_index;
    } else {
      for (unsigned int i = 1; i < this->cols[col_index].p.size(); i += 2) {
        auto len_tmp = 0;
        if (i == this->cols[col_index].p.size() - 1) {
          len_tmp = this->height - this->cols[col_index].p[i];
        } else {
          len_tmp = this->cols[col_index].p[i + 1] - this->cols[col_index].p[i];
        }
        if (len_tmp > best_run_len) {
          best_run_len = len_tmp;
          best_run_index = i;
        }
      }
      return best_run_index;
    }
    return 0;
  }

  /**
   * @brief function to compute length of a run
   * @param col_index column index
   * @param curr_run run index
   */
  unsigned int run_length(unsigned int col_index, unsigned int curr_run) {
    auto len_tmp = 0;
    if (curr_run == this->cols[col_index].p.size() - 1) {
      len_tmp = this->height - this->cols[col_index].p[curr_run];
    } else {
      len_tmp = this->cols[col_index].p[curr_run + 1] -
                this->cols[col_index].p[curr_run];
    }
    return len_tmp;
  }

  /**
   * @brief function to compute k matching statistics with a given query using
   * thresholds
   * @param query haplotype query as std::string
   * @param extend_matches bool to check if extend matching statistics
   * matches with rows
   * @param verbose bool for extra prints
   * @return matching statistics matches
   */
  std::pair<ms_k, std::vector<unsigned int>>
  compute_ms_k(const std::string &query, bool extend_matches = false,
               bool verbose = false) {

    // TODO complete comments
    // compute the match iff |query| is equal to the width of the panel
    if (query.size() != this->width) {
      std::cout << query.size() << " != " << this->width << "\n";
      throw NotEqualLengthException{};
    }

    if (extend_matches && !this->is_extended) {
      // this->extend();
    }

    // initialize matching statistics
    ms_k ms(query.size());
    std::vector<unsigned int> ms_supp(query.size());

    // algorithm begin from the last row of the first column
    // so we obtain the prefix array value (from the samples), the run index
    // and the relative symbol
    unsigned int curr_run = best_run(0, query[0]);
    unsigned curr_pos =
        static_cast<unsigned int>(this->cols[0].sample_end[curr_run]);
    unsigned int curr_index = curr_run != this->cols[0].p.size() - 1
                                  ? this->cols[0].p[curr_run + 1] - 1
                                  : this->height - 1;
    //            std::cout << "curr_pos " << curr_pos << "\n";

    char symbol = get_next_char(this->cols[0].zero_first, curr_run);
    unsigned int s_index = this->cols[0].i_e_k[curr_run];
    bool reset = true;
    // verbose = true;
    //  auto s_index = lf(0, curr_index, query[0]) -
    //  this->cols[0].i_e_k[curr_run];
    //   auto s_index = this->cols[0].p[this->cols[0].sample_end.size() - 1];
    //    iterate over every query's symbol/column index
    for (unsigned int i = 0; i < query.size(); i++) {
      // std::cout << "processed " << i << "\r";
      if (verbose) {
        std::cout << "\nat " << i << " with " << query[i] << " : "
                  << " threshold  " << this->cols[i].t[curr_run] << "\n";
        std::cout << "curr index " << curr_index << " curr run " << curr_run
                  << " curr pos " << curr_pos << " symb " << symbol << "\n";
        std::cout << "sindex " << s_index << " run sindex "
                  << index_to_run(s_index, i) << " with symbol "
                  << get_next_char(this->cols[i].zero_first,
                                   index_to_run(s_index, i))
                  << "\n";
        std::cout << "sindex + k-1 " << s_index + k_smem - 1
                  << " run sindex + k-1 "
                  << index_to_run(s_index + k_smem - 1, i) << "\n";

        std::cout << "zeros " << this->cols[i].count_0 << "\n";
        std::cout << "ones" << this->height - this->cols[i].count_0 << "\n";
      }
      // a lot of cases:
      // - if the pointer in the RLPBWT match the symbol in the query
      //   we proceed simply using lf-mapping (if we are not at the end)
      // - if the pointer in the RLPBWT mismatch the symbol in the query
      //   and we have only that symbol in the column we restart from  the
      //   next column at last index whit the relative prefix array value
      // - otherwise we proceeed using thresholds to select the best
      //   symbol, between the previous and next good symbol (if the
      //   exists), to jump
      //
      //
      auto r_s = index_to_run(s_index, i);
      if (reset == false && query[i] == symbol &&
          ((s_index + k_smem - 1 < this->height &&
            query[i] == get_next_char(this->cols[i].zero_first, r_s) &&
            r_s == index_to_run(s_index + k_smem - 1, i))
           // ||
           // (query[i] == get_next_char(this->cols[i].zero_first,
           //                            index_to_run(s_index, i)) &&
           //  reset == true && run_length(i, curr_run) >= k_smem)
           )) {
        if (verbose) {
          std::cout << "match full:\n";
        }
        ms.row[i] = curr_pos;
        // ms.row_supp[i] = s_index;
        ms_supp[i] = curr_index;

        if (i == 0) {
          ms.len[i] = 1;
          ms.len_supp[i] = 1;
        } else {
          ms.len[i] = ms.len[i - 1] + 1;
          ms.len_supp[i] = ms.len_supp[i - 1] + 1;
        }
        // ms.len[i] = std::min(ms.len[i], ms.len_supp[i]);
        if (i != query.size() - 1) {
          curr_index = lf(i, curr_index, query[i]);
          s_index = lf(i, s_index, query[i]);
          curr_run = index_to_run(curr_index, i + 1);
          symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
          if (verbose) {
            std::cout << "new: " << curr_index << " " << s_index << " "
                      << curr_run << " " << curr_pos << " " << symbol << "\n";
          }
          reset = false;
        }
      } else if ((query[i] == '0' && this->cols[i].count_0 < k_smem) ||
                 (query[i] == '1' &&
                  this->height - this->cols[i].count_0 < k_smem) ||
                 (s_index + k_smem - 1 >= this->height &&
                  query[i] != get_next_char(this->cols[i].zero_first, r_s))) {
        if (verbose) {
          std::cout << "no enough k mismatch/oob k\n";
        }
        // report in matching statistics row vector using panel
        // height as sentinel
        ms.row[i] = this->height;
        // ms.row_supp[i] = this->height;
        ms_supp[i] = this->height;
        ms.len[i] = 0;
        ms.len_supp[i] = 0;
        // update index, run, symbol (as explained before) if we are
        // not at the end
        if (i != query.size() - 1) {
          // std::cout << "start_run"<< std::endl;
          curr_run = best_run(i + 1, query[i + 1]);
          // std::cout << "end_run"<< std::endl;
          // std::cout << std::endl;
          //   std::cout << "start_pos " << i << ""<< std::endl;
          //   std::cout << curr_run << " " << this->cols[i +
          //   1].sample_end.size() << std::endl;
          curr_pos =
              static_cast<unsigned int>(this->cols[i + 1].sample_end[curr_run]);
          //      std::cout << "end_pos"<< std::endl;
          //       std::cout << std::endl;
          //       std::cout << "start_index " << i << "  " << this->width <<
          //       ""<< std::endl;

          curr_index = curr_run != this->cols[i + 1].p.size() - 1
                           ? this->cols[i + 1].p[curr_run + 1] - 1
                           : this->height - 1;
          // std::cout << "end_index"<< std::endl;
          // std::cout << std::endl;
          // std::cout << "start_verbose"<< std::endl;
          if (verbose) {

            std::cout << lf(i + 1, curr_index, query[i + 1]) << " - "
                      << this->cols[i + 1].i_e_k[curr_run] << " = "
                      << lf(i + 1, curr_index, query[i + 1]) -
                             this->cols[i + 1].i_e_k[curr_run]
                      << "\n";
          }
          //     std::cout << "end_verbose"<< std::endl;
          //     std::cout << std::endl;
          //     std::cout << "start_sindex"<< std::endl;
          s_index = lf(i + 1, curr_index, query[i + 1]) -
                    this->cols[i + 1].i_e_k[curr_run];
          // std::cout << "end_sindex"<< std::endl;
          //     std::cout << std::endl;
          // std::cout << "start_symbol"<< std::endl;
          symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
          //  std::cout << "end_symbol"<< std::endl;
          // std::cout << std::endl;
          reset = true;
          if (verbose) {
            std::cout << "update: " << curr_index << " " << curr_pos << " "
                      << symbol << "\n";
            std::cout << "sindex " << s_index << " run sindex "
                      << index_to_run(s_index, i) << "\n";
          }
        }
      } else {
        if (query[i] != symbol) {
          if (verbose) {
            std::cout << "mismatch\n";
          }
          auto thr = this->cols[i].t[curr_run];

          // complete mismatch
          if (this->cols[i].sample_beg.size() == 1 || query[i] == '2') {
            if (verbose) {
              std::cout << "complete mismatch\n";
            }
            // report in matching statistics row vector using panel
            // height as sentinel
            ms.row[i] = this->height;
            // ms.row_supp[i] = this->height;
            ms_supp[i] = this->height;
            ms.len[i] = 0;
            ms.len_supp[i] = 0;
            // update index, run, symbol (as explained before) if we are
            // not at the end
            if (i != query.size() - 1) {
              curr_run = best_run(i + 1, query[i + 1]);
              curr_pos = static_cast<unsigned int>(
                  this->cols[i + 1].sample_end[curr_run]);
              curr_index = curr_run != this->cols[i + 1].p.size() - 1
                               ? this->cols[i + 1].p[curr_run + 1] - 1
                               : this->height - 1;

              symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);

              s_index = lf(i + 1, curr_index, query[i + 1]) -
                        this->cols[i + 1].i_e_k[curr_run];
              reset = true;
              if (verbose) {
                std::cout << "update: " << curr_index << " " << curr_pos << " "
                          << symbol << "\n";
              }
            }
          } else if (curr_run != 0 &&
                     ((curr_index < thr) ||
                      curr_run == this->cols[i].sample_beg.size() - 1)) {
            // if we are above the threshold we go up (if we are not in
            // the first run). We also go up if we are in the last run
            if (verbose) {
              std::cout << "mismatch_up: ";
            }
            curr_index = this->cols[i].p[curr_run] - 1;
            curr_pos = static_cast<unsigned int>(
                this->cols[i].sample_end[curr_run - 1]);
            if (verbose) {
              std::cout << "update: " << curr_index << " " << curr_pos << " "
                        << symbol << "\n";
            }
            // report in matching statistics row vector
            ms.row[i] = curr_pos;
            // ms.row_supp[i] = lf(i, s_index, query[i]);
            if (verbose) {
              std::cout << "curr_index " << curr_index << "i_k "
                        << this->cols[i].i_e_k[index_to_run(curr_index, i)]
                        << "\n";
            }
            // ms.row_supp[i] =
            //    lf(i, curr_index, query[i]) - this->cols[i].i_e_k[curr_run -
            //    1];
            ms_supp[i] = curr_index;
            int tmp_index = (int)i;
            unsigned int len = 0;
            auto rlf = curr_index;

            while (tmp_index >= 0 &&
                   query[tmp_index] ==
                       get_next_char(this->cols[tmp_index].zero_first,
                                     index_to_run(rlf, tmp_index))) {
              if (tmp_index > 0) {
                rlf = reverse_lf(tmp_index, rlf, false);
              }
              tmp_index--;
              len++;
            }

            ms.len[i] = len;
            ms.len_supp[i] = this->cols[i].l_e_k[curr_run - 1];
            // ms.len[i] = std::min(len, ms.len_supp[i]);
            //  update index, run, symbol if we are not at the end
            if (i != query.size() - 1) {
              curr_index = lf(i, curr_index, query[i]);
              if (verbose) {
                std::cout << "update sindex " << curr_index << " "
                          << this->cols[i].i_e_k[curr_run - 1] << "\n";
              }
              s_index = curr_index - this->cols[i].i_e_k[curr_run - 1];
              // s_index = ms_supp[i];
              curr_run = index_to_run(curr_index, i + 1);
              symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
              reset = false;
              if (verbose) {
                std::cout << "new: " << curr_index << " " << curr_run << " "
                          << curr_pos << " " << symbol << "\n";
              }
            }
          } else {
            // we are below threshold  so we go down
            if (verbose) {
              std::cout << "mismatch_down: ";
            }
            curr_index = this->cols[i].p[curr_run + 1];

            curr_pos = static_cast<unsigned int>(
                this->cols[i].sample_beg[curr_run + 1]);

            // report in matching statistics row vector
            ms.row[i] = curr_pos;
            // ms.row_supp[i] = lf(i, curr_index, query[i]) -
            //                 this->cols[i].i_k[index_to_run(curr_index, i)];
            ms_supp[i] = curr_index;
            int tmp_index = (int)i;
            unsigned int len = 0;
            auto rlf = curr_index;
            while (tmp_index >= 0 &&
                   query[tmp_index] ==
                       get_next_char(this->cols[tmp_index].zero_first,
                                     index_to_run(rlf, tmp_index))) {
              if (tmp_index > 0) {
                rlf = reverse_lf(tmp_index, rlf, false);
              }
              tmp_index--;
              len++;
            }
            ms.len[i] = len;
            ms.len_supp[i] = this->cols[i].l_k[curr_run + 1];
            // ms.len[i] = std::min(len, ms.len_supp[i]);
            if (verbose) {
              std::cout << "update: " << curr_index << " " << curr_pos << " "
                        << symbol << "\n";
            }
            // update index, run, symbol if we are not at the end
            if (i != query.size() - 1) {
              curr_index = lf(i, curr_index, query[i]);
              s_index = curr_index - this->cols[i].i_k[curr_run + 1];
              // s_index = ms_supp[i];
              curr_run = index_to_run(curr_index, i + 1);
              symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
              reset = false;
              if (verbose) {
                std::cout << "new: " << curr_index << " " << curr_run << " "
                          << curr_pos << " " << symbol << "\n";
              }
            }
          }
        } else {
          int extra = (s_index + k_smem - 1) - (this->height - 1);
          if (verbose) {
            std::cout << "match not full:\n";
            std::cout << symbol << " vs " << query[i] << "\n";
            std::cout << "curr index " << curr_index << " curr run " << curr_run
                      << " curr pos " << curr_pos << " sizerun "
                      << this->cols[i].p.size() - 1 << "\n";
          }

          if ((curr_index == this->height - 1) ||
              (curr_run != this->cols[i].p.size() - 1 &&
               curr_index == this->cols[i].p[curr_run + 1] - 1)) {
            if (verbose) {
              std::cout << "match end_e\n";
            }
            ms.row[i] = curr_pos;
            // ms.row_supp[i] =
            //    lf(i, curr_index, query[i]) - this->cols[i].i_e_k[curr_run];
            // ms.row_supp[i] = s_index;
            ms_supp[i] = curr_index;

            if (i == 0) {
              ms.len[i] = 1;
              ms.len_supp[i] = 1;
            } else {
              ms.len[i] = ms.len[i - 1] + 1;
              ms.len_supp[i] = this->cols[i].l_e_k[curr_run];
              // ms.len_supp[i] = ms.len_supp[i - 1] + 1;
            }
            // ms.len[i] = std::min(ms.len[i], ms.len_supp[i]);
            if (i != query.size() - 1) {
              curr_index = lf(i, curr_index, query[i]);

              s_index = curr_index - this->cols[i].i_e_k[curr_run];
              curr_run = index_to_run(curr_index, i + 1);
              symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
              reset = false;
              if (verbose) {
                std::cout << "new: " << curr_index << " " << s_index << " "
                          << curr_run << " " << curr_pos << " " << symbol
                          << "\n";
              }
            }

          } else if (extra > 0 && r_s == this->cols[i].p.size() - 1 &&
                     get_next_char(this->cols[i].zero_first, r_s) == query[i] &&
                     r_s == index_to_run(s_index + k_smem - extra - 1, i)) {
            if (r_s <= 2) {
              if (verbose) {
                std::cout << "not enough range no jump\n";
              }
              curr_index = this->height - 1;
              if (curr_run != this->cols[i].p.size() - 1) {
                curr_index = this->cols[i].p[curr_run + 1] - 1;
              }
              ms.row[i] =
                  static_cast<unsigned int>(this->cols[i].sample_end[curr_run]);

              ms_supp[i] = this->cols[i].p[curr_run];

              int tmp_index = (int)i;
              unsigned int len = 0;
              auto rlf = curr_index;
              while (tmp_index >= 0 &&
                     query[tmp_index] ==
                         get_next_char(this->cols[tmp_index].zero_first,
                                       index_to_run(rlf, tmp_index))) {
                if (tmp_index > 0) {
                  rlf = reverse_lf(tmp_index, rlf, false);
                }
                tmp_index--;
                len++;
              }

              ms.len[i] = len;
              ms.len_supp[i] = this->cols[i].l_e_k[curr_run];

              // update index, run, symbol (as explained before) if we are
              // not at the end
              if (i != query.size() - 1) {
                curr_index = lf(i, curr_index, query[i]);
                s_index = curr_index - this->cols[i].i_e_k[curr_run];

                // s_index = ms_supp[i];
                curr_run = index_to_run(curr_index, i + 1);
                symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);

                reset = false;
                if (verbose) {
                  std::cout << "new: " << curr_index << " " << curr_run << " "
                            << curr_pos << " " << symbol << "\n";
                }
                if (verbose) {
                  std::cout << "update: " << curr_index << " " << curr_pos
                            << " " << symbol << "\n";
                }
              }
            } else {
              if (verbose) {
                std::cout << "not enough range yes jump\n";
              }
              curr_run = curr_run - 2;
              curr_index = this->cols[i].p[curr_run + 1] - 1;
              ms.row[i] =
                  static_cast<unsigned int>(this->cols[i].sample_end[curr_run]);
              // ms.row_supp[i] = lf(i, curr_index, query[i]) -
              //                 this->cols[i].i_e_k[index_to_run(curr_index,
              //                 i)];
              ms_supp[i] = this->cols[i].p[curr_run];
              int tmp_index = (int)i;
              unsigned int len = 0;
              auto rlf = curr_index;
              while (tmp_index >= 0 &&
                     query[tmp_index] ==
                         get_next_char(this->cols[tmp_index].zero_first,
                                       index_to_run(rlf, tmp_index))) {
                if (tmp_index > 0) {
                  rlf = reverse_lf(tmp_index, rlf, false);
                }
                tmp_index--;
                len++;
              }
              ms.len[i] = len;
              ms.len_supp[i] = this->cols[i].l_e_k[curr_run];
              // update index, run, symbol (as explained before) if we are
              // not at the end
              if (i != query.size() - 1) {
                curr_index = lf(i, curr_index, query[i]);
                s_index = curr_index - this->cols[i].i_e_k[curr_run + 1];
                // s_index = ms_supp[i];
                curr_run = index_to_run(curr_index, i + 1);
                symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
                reset = false;

                if (verbose) {
                  std::cout << "update: " << curr_index << " " << curr_pos
                            << " " << symbol << "\n";
                }
              }
            }
          } else {
            if (verbose) {
              std::cout << "match not end:\n";
            }
            ms.row[i] = curr_pos;
            // ms.row_supp[i] = s_index;
            ms_supp[i] = curr_index;

            if (i == 0) {
              ms.len[i] = 1;
              ms.len_supp[i] = 1;
            } else {
              ms.len[i] = ms.len[i - 1] + 1;
              ms.len_supp[i] = this->cols[i].l_e_k[curr_run];
            }

            unsigned int b = this->height - 1;
            unsigned int b_i = this->height - 1;
            unsigned int l_i = 0;
            unsigned int r_b = index_to_run(s_index, i);

            unsigned int r_e_b = index_to_run(s_index + this->k_smem - 1, i);

            bool end = false;
            if (verbose) {
              std::cout << "from run " << r_b << " to " << r_e_b << "\n";
            }
            // auto n_r = r_b;
            if (get_next_char(this->cols[i].zero_first, r_b) == query[i]) {
              r_b++;
            }
            if (verbose) {
              std::cout << "update from run " << r_b << " to " << r_e_b << "\n";
            }
            if (r_b != 0 && ((this->cols[i].p[r_b] < this->cols[i].t[r_b]) ||
                             r_b == this->cols[i].sample_beg.size() - 1)) {
              if (verbose) {
                std::cout << "go up\n";
              }
              b_i = this->cols[i].p[r_b] - 1;
              b = this->cols[i].sample_end[r_b - 1];
              l_i = this->cols[i].l_e_k[r_b - 1];
              end = true;
            } else {
              if (verbose) {
                std::cout << "go down\n";
              }
              b_i = this->cols[i].p[r_b + 1];
              b = this->cols[i].sample_beg[r_b + 1];
              l_i = this->cols[i].l_k[r_b + 1];
            }
            if (verbose) {
              std::cout << "b " << b << " b_i " << b_i << " l_i " << l_i
                        << "\n";
            }
            ms.len[i] = std::min(ms.len[i], lce(curr_index, b_i, i, ms.len[i]));
            auto r_b_i = index_to_run(b_i, i);
            ms.row[i] = b;
            curr_pos = b;

            ms.len_supp[i] = l_i;
            ms_supp[i] = b_i;
            if (i != query.size() - 1) {
              curr_index = lf(i, b_i, query[i]);

              if (end) {
                s_index = curr_index - this->cols[i].i_e_k[r_b_i];
              } else {
                s_index = curr_index - this->cols[i].i_k[r_b_i];
              }

              if (verbose) {
                std::cout << "curr " << curr_index << " offset "
                          << this->cols[i].i_e_k[r_b_i] << "\n";

                std::cout << "sindex " << s_index << "\n";
                std::cout << "ri " << r_b_i << "\n";
              }

              curr_run = index_to_run(curr_index, i + 1);
              symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
              reset = false;
              if (verbose) {
                std::cout << "new: " << curr_index << " " << s_index << " "
                          << curr_run << " " << curr_pos << " " << symbol
                          << "\n";
              }
            }
          }
        }
      }

      if (verbose) {
        std::cout << "ms row: " << ms.row[i] << ", ms supp: " << ms_supp[i]
                  << ", ms len: " << ms.len[i]
                  << ", ms len supp: " << ms.len_supp[i]
                  << ", curr index: " << curr_index << ", sindex : " << s_index

                  << ", with size: "
                  << extend_single(ms.row[i],
                                   std::min(ms.len[i], ms.len_supp[i]), i,
                                   ms_supp[i])
                         .size()
                  << ", f_len: " << std::min(ms.len[i], ms.len_supp[i])
                  << " at " << i << " [";
        for (auto x :
             extend_single(ms.row[i], std::min(ms.len[i], ms.len_supp[i]), i,
                           ms_supp[i])) {
          std::cout << x << " ";
        }
        std::cout << "\n";
      }
    }
    // exit(0);
    return std::make_pair(ms, ms_supp);
  }

  /**
   * @brief function to compute matching statistics with a given query using
   * thresholds
   * @param query haplotype query as std::string
   * @param extend_matches bool to check if extend matching statistics matches
   * with rows
   * @param verbose bool for extra prints
   * @return matching statistics matches
   */
  std::pair<ms, std::vector<unsigned int>>
  compute_ms(const std::string &query, bool extend_matches = false,
             bool verbose = false) {

    // compute the match iff |query| is equal to the width of the panel
    if (query.size() != this->width) {
      std::cout << query.size() << " != " << this->width << "\n";
      throw NotEqualLengthException{};
    }
    // if required extend with the phi support struct (iff not already
    // extended)
    if (extend_matches && !this->is_extended) {
      // this->extend();
    }
    // initialize matching statistics
    ms ms(query.size());
    std::vector<unsigned int> ms_supp(query.size());

    // algorithm begin from the last row of the first column
    // so we obtain the prefix array value (from the samples), the run index
    // and the relative symbol
    auto curr_pos = static_cast<unsigned int>(
        this->cols[0].sample_end[this->cols[0].sample_end.size() - 1]);
    auto curr_index = curr_pos;
    //            std::cout << "curr_pos " << curr_pos << "\n";
    unsigned int curr_run = index_to_run(curr_index, 0);
    char symbol = get_next_char(this->cols[0].zero_first, curr_run);
    // iterate over every query's symbol/column index
    for (unsigned int i = 0; i < query.size(); i++) {

      // std::cout << "processed " << i << "\r";
      if (verbose) {
        std::cout << "at " << i << ": " << curr_run << " "
                  << this->cols[i].t[curr_run] << "\n";
        std::cout << curr_index << " " << curr_run << " " << curr_pos << " "
                  << symbol << "\n";
      }
      // a lot of cases:
      // - if the pointer in the RLPBWT match the symbol in the query
      //   we proceed simply using lf-mapping (if we are not at the end)
      // - if the pointer in the RLPBWT mismatch the symbol in the query
      //   and we have only that symbol in the column we restart from  the
      //   next column at last index whit the relative prefix array value
      // - otherwise we proceeed using thresholds to select the best
      //   symbol, between the previous and next good symbol (if the
      //   exists), to jump
      if (query[i] == symbol) {
        if (verbose) {
          std::cout << "match:\n";
        }
        // report in matching statistics row vector
        ms.row[i] = curr_pos;
        ms_supp[i] = curr_index;
        // update index, run, symbol if we are not at the end
        if (i != query.size() - 1) {
          curr_index = lf(i, curr_index, query[i]);

          curr_run = index_to_run(curr_index, i + 1);
          symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
          if (verbose) {
            std::cout << "new: " << curr_index << " " << curr_run << " "
                      << curr_pos << " " << symbol << "\n";
          }
        }
      } else {
        // get threshold
        auto thr = this->cols[i].t[curr_run];

        if (this->cols[i].sample_beg.size() == 1 || query[i] == '2') {
          if (verbose) {
            std::cout << "complete mismatch\n";
          }
          // report in matching statistics row vector using panel
          // height as sentinel
          ms.row[i] = this->height;
          ms_supp[i] = this->height;

          // update index, run, symbol (as explained before) if we are
          // not at the end
          if (i != query.size() - 1) {
            curr_pos = static_cast<unsigned int>(
                this->cols[i + 1]
                    .sample_end[this->cols[i + 1].sample_end.size() - 1]);
            curr_index = this->height - 1;

            curr_run = index_to_run(curr_index, i + 1);
            symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
            if (verbose) {
              std::cout << "update: " << curr_index << " " << curr_pos << " "
                        << symbol << "\n";
            }
          }
        } else if (curr_run != 0 &&
                   ((curr_index < thr) ||
                    curr_run == this->cols[i].sample_beg.size() - 1)) {
          // if we are above the threshold we go up (if we are not in
          // the first run). We also go up if we are in the last run
          if (verbose) {
            std::cout << "mismatch_up: ";
          }
          curr_index = this->cols[i].p[curr_run] - 1;
          curr_pos =
              static_cast<unsigned int>(this->cols[i].sample_end[curr_run - 1]);
          if (verbose) {
            std::cout << "update: " << curr_index << " " << curr_pos << " "
                      << symbol << "\n";
          }
          // report in matching statistics row vector
          ms.row[i] = curr_pos;

          ms_supp[i] = curr_index;

          // update index, run, symbol if we are not at the end
          if (i != query.size() - 1) {
            curr_index = lf(i, curr_index, query[i]);
            curr_run = index_to_run(curr_index, i + 1);
            symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
            if (verbose) {
              std::cout << "new: " << curr_index << " " << curr_run << " "
                        << curr_pos << " " << symbol << "\n";
            }
          }
        } else {
          // we are below threshold  so we go down
          if (verbose) {
            std::cout << "mismatch_down: ";
          }
          curr_index = this->cols[i].p[curr_run + 1];

          curr_pos =
              static_cast<unsigned int>(this->cols[i].sample_beg[curr_run + 1]);

          // report in matching statistics row vector
          ms.row[i] = curr_pos;
          ms_supp[i] = curr_index;

          if (verbose) {
            std::cout << "update: " << curr_index << " " << curr_pos << " "
                      << symbol << "\n";
          }
          // update index, run, symbol if we are not at the end
          if (i != query.size() - 1) {
            curr_index = lf(i, curr_index, query[i]);
            curr_run = index_to_run(curr_index, i + 1);
            symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
            if (verbose) {
              std::cout << "new: " << curr_index << " " << curr_run << " "
                        << curr_pos << " " << symbol << "\n";
            }
          }
        }
      }
    }

    // compute the len vector using random access on the panel, proceeding
    // from left to right

    for (unsigned int i = 0; i < ms.len.size(); i++) {
      if (ms.row[i] == this->height) {
        // if we have the sentinel in row vector than the length is 0
        ms.len[i] = 0;
      } else if (i != 0 && ms.row[i] == ms.row[i - 1] && ms.len[i - 1] != 0) {
        ms.len[i] = ms.len[i - 1] + 1;
      } else {

        int tmp_index = (int)i;
        unsigned int len = 0;
        auto rlf = ms_supp[tmp_index];
        while (tmp_index >= 0 &&
               query[tmp_index] ==
                   get_next_char(this->cols[tmp_index].zero_first,
                                 index_to_run(rlf, tmp_index))) {
          if (tmp_index > 0) {
            rlf = reverse_lf(tmp_index, rlf, false);
          }
          tmp_index--;
          len++;
        }
        ms.len[i] = len;
      }
    }
    if (verbose) {
      std::cout << ms << "\n";
    }

    return std::make_pair(ms, ms_supp);
  }

  std::vector<unsigned int> intersection(std::vector<unsigned int> v1,
                                         std::vector<unsigned int> v2) {
    std::vector<unsigned int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                          back_inserter(v3));
    return v3;
  }

  std::vector<std::pair<unsigned int, unsigned int>>
  intersection_pairs(std::vector<std::pair<unsigned int, unsigned int>> v1,
                     std::vector<std::pair<unsigned int, unsigned int>> v2) {
    std::vector<std::pair<unsigned int, unsigned int>> v3;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                          back_inserter(v3));
    return v3;
  }

  std::vector<std::pair<unsigned int, unsigned int>>
  filter_comp(std::vector<std::pair<unsigned int, unsigned int>> p,
              std::vector<unsigned int> r) {
    std::vector<std::pair<unsigned int, unsigned int>> res;
    std::set<unsigned int> s(r.begin(), r.end());
    for (auto e : p) {
      if (s.find(e.first) != s.end() && s.find(e.second) != s.end())
        res.push_back({e});
    }
    return res;
  }
  std::vector<unsigned int>
  red_pairs(std::vector<std::pair<unsigned int, unsigned int>> p) {
    std::vector<unsigned int> res;
    for (auto pp : p) {
      res.push_back(pp.first);
      res.push_back(pp.second);
    }
    std::unordered_set<unsigned int> s;
    for (int i : res)
      s.insert(i);
    res.assign(s.begin(), s.end());
    sort(res.begin(), res.end());
    return res;
  }

  std::vector<std::pair<unsigned int, unsigned int>>
  generate_pairs(const std::vector<unsigned int> &input) {
    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (size_t i = 0; i < input.size(); ++i) {
      for (size_t j = i + 1; j < input.size(); ++j) {
        if (input[i] < input[j])
          result.emplace_back(input[i], input[j]);
        else
          result.emplace_back(input[j], input[i]);
      }
    }

    return result;
  }

  void
  query_match(const char *ref, const char *filename, const char *out,
              std::unordered_set<unsigned int> ps,
              std::unordered_set<unsigned int> ps2,
              std::vector<std::string> samples_t,
              std::vector<std::string> samples, std::vector<unsigned int> sites,
              std::vector<std::vector<std::string>> data_t,
              std::vector<std::vector<std::string>> data,
              std::tuple<std::vector<sdsl::sd_vector<>>,
                         std::vector<std::string>, std::vector<unsigned int>>
                  ref_data,
              bool unsafe, bool v = false) {
    std::ofstream out_match(out);
    std::string new_row;
    std::string line;
    std::vector<std::string> queries_panel;
    auto ref_panel = std::get<0>(ref_data);
    auto ref_samples = std::get<1>(ref_data);
    auto ref_sites = std::get<2>(ref_data);
    // while (getline(input_matrix, line) && !line.empty()) {
    //   std::istringstream is_col(line);
    //   is_col >> new_row;
    //   queries_panel.push_back(new_row);
    // }

    htsFile *fp = hts_open(filename, "rb");
    if (fp == NULL) {
      throw FileNotFoundException{};
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec = bcf_init();

    unsigned int count = 0;
    std::string last_col;
    std::string last_column;
    std::string new_column;
    std::string tmp_column;
    while (bcf_read(fp, hdr, rec) >= 0) {
      bcf_unpack(rec, BCF_UN_ALL);
      // if (auto search = ps.find(rec->pos); search == ps.end())
      //   continue;
      new_column = "";
      int32_t *gt_arr = NULL, ngt_arr = 0;
      int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
      ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
      int max_ploidy = ngt / nsmpl;
      for (i = 0; i < nsmpl; i++) {
        int32_t *ptr = gt_arr + i * max_ploidy;
        for (j = 0; j < max_ploidy; j++) {
          if (ptr[j] == bcf_int32_vector_end)
            break;

          if (bcf_gt_is_missing(ptr[j]))
            exit(-1);

          int allele_index = bcf_gt_allele(ptr[j]);
          if (j == 0) {
            new_column += std::to_string(allele_index);
          } else if (j == 1 &&
                     new_column[new_column.size() - 1] - '0' != allele_index) {
            new_column[new_column.size() - 1] = '2';
          }
        }
      }

      queries_panel.push_back(new_column);
      free(gt_arr);
    }
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    bcf_destroy(rec);
    std::string query;
    std::vector<std::string> queries;
    for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
      for (auto &j : queries_panel) {
        query.push_back(j[i]);
      }
      queries.push_back(query);
      query.clear();
    }

    auto n_queries = queries.size();

    std::vector<ms_matches> matches_vec(n_queries);
    std::vector<std::vector<unsigned int>> phased_haplos;
    std::vector<std::vector<bool>> is_phases_full;
    // bool unsafe = false;
    for (unsigned int i = 0; i < n_queries; i++) {
      // auto x = this->match_thr(queries[i]);
      // std::cout << x;
      matches_vec[i] = this->left_mpsc(queries[i], true);
      // std::cout << matches_vec[i];
      unsigned int j = matches_vec[i].basic_matches.size() - 1;
      std::vector<unsigned int> r_rows;
      std::vector<unsigned int> l_rows;
      unsigned int l_col;
      unsigned int r_col;
      std::vector<unsigned int> rows;
      std::vector<std::vector<unsigned int>> mat;
      std::vector<std::pair<unsigned int, unsigned int>> comp;
      std::vector<std::pair<unsigned int, unsigned int>> comp_supp;
      std::vector<std::pair<std::vector<std::pair<unsigned int, unsigned int>>,
                            unsigned int>>
          recomb;
      bool is_recomb = false;
      unsigned int len_pref = 0;
      unsigned int len_suff = 0;
      int k = 0;
      while (k < queries[i].size() && queries[i][k] == '2') {
        len_pref++;
        k++;
      }
      k = queries[i].size() - 1;

      while (k >= 0 && queries[i][k] == '2') {
        len_suff++;
        k--;
      }

      // r_rows = matches_vec[i].haplos[matches_vec[i].haplos.size() - 1];
      // l_rows = matches_vec[i].haplos[matches_vec[i].haplos.size() - 2];
      // rows = intersection(r_rows, l_rows);
      // l_col =
      // std::get<2>(matches_vec[i].basic_matches[matches_vec[i].basic_matches.size()
      // - 1]) + 1; r_col =
      // std::get<2>(matches_vec[i].basic_matches[matches_vec[i].basic_matches.size()
      // - 2]) -
      // std::get<1>(matches_vec[i].basic_matches[matches_vec[i].basic_matches.size()
      // - 2]); mat = get_sub_matrix(rows, l_col, r_col); comp =
      // get_comp(mat); l_rows = red_pairs(comp);
      //
      std::vector<bool> is_safety_phased(this->width, true);
      std::vector<std::pair<unsigned int, unsigned int>> unsafe_pos;

      l_rows = matches_vec[i].haplos[j];
      if (len_pref != 0) {
        std::vector<unsigned int> v(this->height);
        std::iota(std::begin(v), std::end(v), 0);
        mat = get_sub_matrix(v, 0, len_pref - 1);
        comp = get_comp(mat, v);

        l_rows = red_pairs(comp);
        unsafe_pos.push_back({0, len_pref - 1});
      }
      // TODO aggiungere switch tra mpsc consecutivi senxa 2 in mezzo
      while (j >= 1) {

        r_rows = matches_vec[i].haplos[j - 1];
        rows = intersection(r_rows, l_rows);
        l_col = std::get<2>(matches_vec[i].basic_matches[j]) + 1;
        r_col = std::get<2>(matches_vec[i].basic_matches[j - 1]) -
                std::get<1>(matches_vec[i].basic_matches[j - 1]);
        // if (j >= matches_vec[i].basic_matches.size() - 100) {
        //   for (auto c : comp)
        //     std::cout << c.first << "," << c.second << " ";
        //   std::cout << std::endl;
        //   std::cout << j << " " << l_col << " " << r_col << std::endl;
        //   for (auto c : l_rows)
        //     std::cout << c << " ";
        //   std::cout << std::endl;
        //   for (auto c : r_rows)
        //     std::cout << c << " ";
        //   std::cout << std::endl;
        //   for (auto c : rows)
        //     std::cout << c << " ";
        //   std::cout << std::endl;
        // }
        if (r_rows.size() == 1) {
          unsafe_pos.push_back({l_col, r_col});
          comp.clear();
          j -= 2;
          if (j <= 0)
            break;
          l_rows = matches_vec[i].haplos[j];
          continue;
        }
        if ((l_col - 1) == r_col) {
          is_recomb = true;
          if (rows.size() >= 2) {
            comp = filter_comp(comp, r_rows);
            l_rows = red_pairs(comp);
            if (comp.empty()) {
              l_rows = matches_vec[i].haplos[j];
              comp = generate_pairs(l_rows);
            }
            if (!comp.empty()) {
              if (recomb.empty()) {
                recomb.push_back(
                    {comp, std::get<2>(matches_vec[i].basic_matches[j - 1])});
              } else {
                recomb[recomb.size() - 1] = {
                    comp, std::get<2>(matches_vec[i].basic_matches[j - 1])};
              }
            } else {
              l_rows = matches_vec[i].haplos[j];
            }
          } else {
            if (l_rows.size() < 2)
              l_rows = matches_vec[i].haplos[j];
            comp = generate_pairs(l_rows);
            if (!comp.empty()) {
              if (recomb.empty()) {
                recomb.push_back(
                    {comp, std::get<2>(matches_vec[i].basic_matches[j - 1])});
              } else {
                recomb[recomb.size() - 1] = {
                    comp, std::get<2>(matches_vec[i].basic_matches[j - 1])};
              }
            } else {
              l_rows = matches_vec[i].haplos[j];
            }
          }
        } else {
          mat = get_sub_matrix(rows, l_col, r_col);
          // TODO se mat vuota intersezione tra pairs e nuove righe (dei pairs
          // tengo solo quelle con righe negli haplosm di MPSC successiva)
          comp_supp = get_comp(mat, rows);

          auto back = comp;
          if ((j == matches_vec[i].basic_matches.size() - 1 && len_pref == 0) ||
              comp.empty()) {
            l_rows = red_pairs(comp_supp);
            comp = comp_supp;
          } else {
            comp = intersection_pairs(comp, comp_supp);
            l_rows = red_pairs(comp);
          }
          // for (auto c : comp)
          //   std::cout << c.first << "," << c.second << " ";
          // std::cout << j << std::endl;
          if (!comp.empty()) {
            if (recomb.empty()) {
              recomb.push_back(
                  {comp, std::get<2>(matches_vec[i].basic_matches[j - 1])});

            } else {
              recomb[recomb.size() - 1] = {
                  comp, std::get<2>(matches_vec[i].basic_matches[j - 1])};
            }
          } else {
            // if (!back.empty() && std::get<2>(matches_vec[i].basic_matches[j])
            // <
            //                          recomb[recomb.size() - 1].second)
            //   recomb.push_back(
            //       {back, std::get<2>(matches_vec[i].basic_matches[j])});
            is_recomb = true;
            std::vector<unsigned int> v(this->height);
            std::iota(std::begin(v), std::end(v), 0);
            mat = get_sub_matrix(v, l_col, r_col);
            comp = get_comp(mat, v);
            l_rows = red_pairs(comp);
            // if (comp.empty()) {

            //   l_rows = matches_vec[i].haplos[j];
            //   comp = generate_pairs(l_rows);
            // }
            if (!unsafe) {
              if (!comp.empty())
                recomb.push_back(
                    {comp, std::get<2>(matches_vec[i].basic_matches[j - 1])});
              else
                l_rows = matches_vec[i].haplos[j];
            }
            if (unsafe) {
              unsafe_pos.push_back({l_col, r_col});
              comp.clear();
              l_rows = matches_vec[i].haplos[j];
              comp = generate_pairs(l_rows);
              if (!comp.empty())
                recomb.push_back(
                    {comp, std::get<2>(matches_vec[i].basic_matches[j - 1])});
            }
          }
        }
        j--;
      }

      if (len_suff != 0) {
        mat = get_sub_matrix(l_rows, this->width - len_suff, this->width - 1);
        comp_supp = get_comp(mat, l_rows);
        if (!unsafe) {
          comp = intersection_pairs(comp, comp_supp);
          if (!comp.empty())
            recomb.push_back({comp, this->width - 1});
          else {
            // comp = generate_pairs(l_rows);
            std::vector<unsigned int> v(this->height);
            std::iota(std::begin(v), std::end(v), 0);
            mat = get_sub_matrix(v, this->width - len_suff, this->width - 1);
            comp = get_comp(mat, v);
            if (!comp.empty())
              recomb.push_back({comp, this->width - 1});
          }
        }
        if (unsafe)
          unsafe_pos.push_back({this->width - len_suff + 1, this->width - 1});
      }
      // std::cout << "computed haps for query: " << i
      //           << "having recomb: " << is_recomb << "\n";
      // if (is_recomb)
      //   std::cout << "Is recomb" << i << std::endl;
      // else
      //   std::cout << "Is not recomb" << i << std::endl;
      if (unsafe) {
        for (auto upos : unsafe_pos) {
          for (unsigned int u = upos.first; u <= upos.second; u++)
            is_safety_phased[u] = false;
        }
      }
      // is_phases_full.push_back(is_safety_phased);
      // std::vector<std::pair<std::vector<std::pair<unsigned int, unsigned
      // int>>,
      //                    unsigned int>>
      // recomb_f;
      // for (auto r_comp : recomb) {

      //   std::cout << "[ ";
      //   for (auto c : r_comp.first) {
      //     std::cout << "(" << c.first << ", " << c.second << ") ";
      //   }
      //   std::cout << ", " << r_comp.second << "]; \n";
      // }
      std::cout << std::endl;
      optimize_recomb(recomb);
      std::vector<unsigned int> target_pos;
      for (auto d : data_t)
        target_pos.push_back(std::stoi(d[1]) - 1);
      // std::cout << comp.size() << " " << recomb.size() << "\n";
      // for (auto c : comp) {
      //   std::cout << "(" << c.first << ", " << c.second << ")\n";
      // }
      // std::cout << std::endl;

      // for (auto r_comp : recomb) {
      //   std::cout << "[ ";
      //   for (auto c : r_comp.first) {
      //     std::cout << "(" << c.first << ", " << c.second << ") ";
      //   }
      //   std::cout << ", " << r_comp.second << "]; \n";
      // }

      // std::cout << std::endl;

      // if(is_recomb && recomb.size() > 1) {
      //      std::cout << "@" << i << std::endl;
      // for (auto r_comp : recomb) {
      //   std::cout << "[ ";
      //   for (auto c : r_comp.first) {
      //     std::cout << "(" << c.first << ", " << c.second << ") ";
      //   }
      //   std::cout << ", \n" << r_comp.second << "]; \n";
      // }

      /* auto f = recomb[0];
       for(unsigned int ii = 1; ii < recomb.size(); ii++) {
         auto tmp = intersection_pairs(f.first, recomb[ii].first);
         if(!tmp.empty() && !recomb_f.empty()){
           recomb_f[recomb_f.size() - 1]={tmp, recomb[ii].second};
         }else{
           recomb_f.push_back(recomb[ii]);
         }
         f = recomb_f[recomb_f.size() - 1];
       }
     }
     recomb = recomb_f;
   */
      // htsFile *fpr2 = bcf_open(ref, "r");

      // bcf_hdr_t *hdrr2 = bcf_hdr_read(fpr2);
      if (!is_recomb) {
        // std::cout << "@" << i << std::endl;
        //  for (auto c : comp) {
        //    std::cout << "(" << c.first << ", " << c.second << ")\n";
        //  }
        // int hap1 = comp[0].first % 2 == 0 ? 0 : 1;
        // auto h1 = gh(ref, samples[int(std::floor(comp[0].first / 2))], hap1,
        // ps,
        //              fpr2, hdrr2);
        // int hap2 = comp[0].second % 2 == 0 ? 0 : 1;
        // auto h2 = gh(ref, samples[int(std::floor(comp[0].second / 2))], hap2,
        //              ps, fpr2, hdrr2);

        // auto h1 = slice_sd_vector(this->panel[comp[0].first], 0,
        // this->width); auto h2 = slice_sd_vector(this->panel[comp[0].second],
        // 0, this->width); phased_haplos.push_back(h1);
        // phased_haplos.push_back(h2);
        std::vector<unsigned int> samples_ref1(this->width, comp[0].first);
        std::vector<unsigned int> samples_ref2(this->width, comp[0].second);
        auto map_pos1 = get_ref_pos(samples_ref1, target_pos, ref_sites);
        auto map_pos2 = get_ref_pos(samples_ref2, target_pos, ref_sites);

        std::vector<unsigned int> tmp1 = getfromref(ref_panel, map_pos1);
        std::vector<unsigned int> tmp2 = getfromref(ref_panel, map_pos2);
        phased_haplos.push_back(tmp1);
        phased_haplos.push_back(tmp2);
        std::vector<bool> t(data.size(), true);
        is_phases_full.push_back(t);
      } else {
        // TODO check existence of recomb
        // std::cout << "@" << i << std::endl;
        std::vector<unsigned int> tmp1;
        std::vector<unsigned int> tmp2;
        unsigned int start = 0;
        std::vector<unsigned int> samples_ref1(this->width);
        std::vector<unsigned int> samples_ref2(this->width);
        // std::cout << recomb.size() << std::endl;
        for (auto r_comp : recomb) {
          // std::cout << start << " " << r_comp.second << " " << this->width
          //           << std::endl;
          if (r_comp.second == this->width)
            r_comp.second--;
          if (start > r_comp.second)
            break;
          // std::cout << start << " " << r_comp.second << " " << this->width
          //           << " " << r_comp.first.size() << std::endl;
          //  int hap1 = r_comp.first[0].first % 2 == 0 ? 0 : 1;

          // auto h1 = gh(ref,
          // samples[int(std::floor(r_comp.first[0].first / 2))],
          //              hap1, ps, fpr2, hdrr2);
          // int hap2 = r_comp.first[0].second % 2 == 0 ? 0 : 1;
          // auto h2 =
          //     gh(ref, samples[int(std::floor(r_comp.first[0].second /
          //     2))],
          //        hap2, ps, fpr2, hdrr2);

          // auto h1 = slice_sd_vector(this->panel[r_comp.first[0].first],
          // start,
          //                           r_comp.second);
          // auto h2 =
          // slice_sd_vector(this->panel[r_comp.first[0].second], start,
          //                           r_comp.second);
          // //
          // tmp1.insert(tmp1.end(), h1.begin(), h1.end());
          // tmp2.insert(tmp2.end(), h2.begin(), h2.end());
          for (int ii = start; ii <= r_comp.second; ii++) {
            // std::cout << r_comp.first[0].first << std::endl;
            // std::cout << r_comp.first[0].second << std::endl;

            // std::cout << "-------" << std::endl;
            samples_ref1[ii] = r_comp.first[0].first;
            // samples[int(std::floor(r_comp.first[0].first / 2))];
            samples_ref2[ii] = r_comp.first[0].second;
            // samples[int(std::floor(r_comp.first[0].second / 2))];
            //  samples_ref1[ii] =
            //      samples[int(std::floor(r_comp.first[0].first / 2))];
            //  samples_ref2[ii] =
            //      samples[int(std::floor(r_comp.first[0].second / 2))];
          }
          // std::cout << "------aaaaaaaa-" << std::endl;
          start = r_comp.second + 1;
          // std::cout << "[ ";
          // for (auto c : r_comp.first) {
          //   std::cout << "(" << c.first << ", " << c.second << ") ";
          // }
          // std::cout << ", " << r_comp.second << "];\n";
        }
        // std::cout << "------nnbbbbbbbba-" << std::endl;
        // std::cout << samples_ref1.size() << " " << target_pos.size() << " "
        //          << ref_sites.size() << "\n";
        // for (auto x : samples_ref1) {
        //  std::cout << x << " ";
        //}
        // std::cout << std::endl;
        auto map_pos1 = get_ref_pos(samples_ref1, target_pos, ref_sites);
        // for (auto m : map_pos1) {
        //   std::cout << std::get<0>(m) << " " << std::get<1>(m) << " "
        //             << std::get<1>(m) << "\n";
        // }
        auto map_pos2 = get_ref_pos(samples_ref2, target_pos, ref_sites);
        tmp1 = getfromref(ref_panel, map_pos1);
        // std::cout << tmp1.size() << "\n";
        tmp2 = getfromref(ref_panel, map_pos2);
        phased_haplos.push_back(tmp1);
        phased_haplos.push_back(tmp2);
        // std::cout << is_safety_phased.size() << " " << target_pos.size() << "
        // "
        //           << ref_sites.size() << "\n";
        is_phases_full.push_back(
            expand_phased(is_safety_phased, target_pos, ref_sites));

        // for (auto x : ref_panel) {
        //   std::cout << x << "\n";
        // }
        // std::cout << std::endl;
      }
      // bcf_hdr_destroy(hdrr2);
      // bcf_close(fpr2);
    }

    // for (auto s : phased_haplos) {
    //   for (auto hh : s) {
    //     std::cout << hh << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << phased_haplos[0].size() << " " << data.size() << "\n";
    auto hs = cp(phased_haplos);

    // for (auto hh : hs) {
    //   for (auto hhh : hh) {
    //     std::cout << "(" << hhh.first << "," << hhh.second << ") ";
    //   }
    //   std::
    //    cout << std::endl;
    // }
    // unsafe = true;
    // std::cout << is_phases_full.size() << " " << hs[0].size() << " \n";
    std::vector<bool> ph(data.size());
    if (unsafe) {
      if (is_phases_full.size() == 1) {
        ph = is_phases_full[0];
      } else {
        for (const auto &v : is_phases_full) {
          for (size_t i = 0; i < is_phases_full[0].size(); ++i) {
            ph[i] = ph[i] && v[i];
          }
        }
      }
    } else {
      for (unsigned int u = 0; u < data.size(); u++) {
        ph[u] = true;
      }
    }
    // std::cout << ph.size() << " " << hs[0].size() << " \n";
    //   for (unsigned int u = 0; u < this->width; u++) {
    //     std::cout << ph[u];
    //   }
    htsFile *fpo = bcf_open(out, "wb");

    bcf_hdr_t *hdro = bcf_hdr_init("w");
    if (!hdro) {
      fprintf(stderr, "Header creation error\n");
      return;
    }
    std::string chra = "";
    std::set<std::string> chrs;
    for (auto d : data) {
      auto tc = d[0];
      if (auto search = chrs.find(tc); search == chrs.end()) {
        chra += (tc + ",");
        chrs.insert(tc);
      }
    }
    chra.erase(chra.size() - 1);
    std::string header = "##contig=<ID=" + chra + ",length=" +
                         std::to_string(ref_sites[ref_sites.size() - 1]) + ">";
    bcf_hdr_append(hdro, "##fileformat=VCFv4.2");
    bcf_hdr_append(hdro, header.c_str());
    bcf_hdr_append(
        hdro,
        "##FORMAT=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">");

    bcf_hdr_append(
        hdro, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    for (auto sam : samples_t) {
      bcf_hdr_add_sample(hdro, sam.c_str());
    }

    bcf_hdr_add_sample(hdro, NULL);

    if (bcf_hdr_write(fpo, hdro) < 0) {
      fprintf(stderr, "Header writing errorr\n");
      return;
    }
    int cl = 0;
    for (int i = 0; i < data.size(); i++) {
      // std::cout << i << "\n";
      bcf1_t *reco = bcf_init();
      reco->rid = bcf_hdr_name2id(hdro, data[i][0].c_str());
      // std::cout << reco->rid << " " << data[i][0].c_str() << "\n";
      reco->pos = std::stoi(data[i][1]) - 1;

      bcf_update_id(hdro, reco, data[i][2].c_str());
      // std::cout << reco->pos << "\n";
      auto al = data[i][3] + ',' + data[i][4];
      // std::cout << al << "\n";
      bcf_update_alleles_str(hdro, reco, al.c_str());

      int32_t gt_arr[samples_t.size() * 2];
      // if (ph[i]) {
      //   // std::cout << "here\n";
      //   for (int j = 0; j < samples_t.size(); j++) {
      //     // std::cout << "writing: " << hs[j][i].first << " " <<
      //     // hs[j][i].second
      //     //           << "\n";
      //     gt_arr[j * 2] = bcf_gt_phased(hs[j][i].first);
      //     gt_arr[j * 2 + 1] = bcf_gt_phased(hs[j][i].second);
      //   }
      // } else {
      //   for (int j = 0; j < samples_t.size(); j++) {
      //     // std::cout << "writing: " << 2 << " at " << i << "\n";
      //     gt_arr[j * 2] = bcf_gt_unphased(0);
      //     gt_arr[j * 2 + 1] = bcf_gt_unphased(1);
      //   }
      // }
      //
      //
      if (ph[i]) {
        int32_t f = 1;
        bcf_update_format_int32(hdro, reco, "IMP", &f, 0);
      }
      if (!ph[i]) {
        cl++;
      }
      // else {
      //   bcf_update_format_int32(hdro, reco, "IMP", 1, 1);
      // }
      for (int j = 0; j < samples_t.size(); j++) {
        gt_arr[j * 2] = bcf_gt_phased(hs[j][i].first);
        gt_arr[j * 2 + 1] = bcf_gt_phased(hs[j][i].second);
      }
      bcf_update_genotypes(hdro, reco, gt_arr, samples_t.size() * 2);

      if (bcf_write(fpo, hdro, reco) < 0) {
        fprintf(stderr, "Writing variant error %d\n", i);
        return;
      }

      bcf_destroy(reco);
    }
    bcf_hdr_destroy(hdro);
    bcf_close(fpo);
    std::cout << "Low prob over target " << cl << std::endl;
  }

  void optimize_recomb(
      std::vector<std::pair<std::vector<std::pair<unsigned int, unsigned int>>,
                            unsigned int>> &recomb) {
    std::map<std::pair<int, int>, int> counts;
    std::set<std::pair<int, int>> used;

    for (const auto &r : recomb) {
      for (const auto &p : r.first) {
        counts[p]++;
      }
    }

    for (size_t i = 0; i < recomb.size(); i++) {
      std::pair<int, int> best = recomb[i].first[0];
      int m = 0;

      for (const auto &p : recomb[i].first) {
        if (counts[p] > m || (counts[p] == m && used.count(p))) {
          best = p;
          m = counts[p];
        }
      }
      recomb[i].first = {best};
      used.insert(best);
    }
  }

  std::vector<unsigned int> getfromref(
      const std::vector<sdsl::sd_vector<>> &matrix,
      const std::vector<std::tuple<unsigned int, size_t, size_t>> &intervals) {

    std::vector<unsigned int> result;

    if (matrix.empty())
      return result;

    unsigned int rowIdx = std::get<0>(intervals[0]);
    size_t startIdx = std::get<1>(intervals[0]);
    size_t endIdx = std::get<2>(intervals[0]);

    for (auto &i : intervals) {
      auto h = slice_sd_vector(matrix[std::get<0>(i)], std::get<1>(i),
                               std::get<2>(i));
      result.insert(result.end(), h.begin(), h.end());
    }

    return result;
  }

  // std::vector<std::tuple<bool, size_t, size_t>>
  // expand_phased(const std::vector<bool> &phased,
  //               const std::vector<unsigned int> &target_pos,
  //               const std::vector<unsigned int> &ref_pos) {

  //   // std::vector<unsigned int> res(ref_pos.size());
  //   // size_t i = 0;
  //   // unsigned int lastValue = 0;

  //   // for (size_t j = 0; j < ref_pos.size(); ++j) {
  //   //   if (i < target_pos.size() && ref_pos[j] == target_pos[i]) {
  //   //     lastValue = samples[i];
  //   //     ++i;
  //   //   }
  //   //   res[j] = lastValue;
  //   // }

  //   // return res;
  //   //
  //   std::vector<std::tuple<bool, size_t, size_t>> intervals;
  //   size_t i = 0;
  //   size_t start_idx = 0;
  //   bool lastValue = phased[0];

  //   for (size_t j = 0; j < ref_pos.size(); ++j) {
  //     if (i + 1 < target_pos.size() && ref_pos[j] >= target_pos[i + 1]) {
  //       intervals.emplace_back(lastValue, start_idx, j);
  //       start_idx = j + 1;
  //       lastValue = phased[++i];
  //     }
  //   }

  //   intervals.emplace_back(lastValue, start_idx, ref_pos.size() - 1);

  //   return compactIntervals(intervals);
  // }
  //
  std::vector<bool> expand_phased(const std::vector<bool> &phased,
                                  const std::vector<unsigned int> &target_pos,
                                  const std::vector<unsigned int> &ref_pos) {

    std::vector<bool> res(ref_pos.size(), true);
    size_t i = 0;
    size_t start_idx = 0;
    bool lastValue = phased[0];

    for (size_t j = 0; j < ref_pos.size(); ++j) {
      if (i + 1 < target_pos.size() && ref_pos[j] == target_pos[i]) {
        res[j] = phased[i];
        i++;
      }
    }

    return res;
  }

  std::vector<std::tuple<unsigned int, size_t, size_t>>
  get_ref_pos(const std::vector<unsigned int> &samples,
              const std::vector<unsigned int> &target_pos,
              const std::vector<unsigned int> &ref_pos) {

    // std::vector<unsigned int> res(ref_pos.size());
    // size_t i = 0;
    // unsigned int lastValue = 0;

    // for (size_t j = 0; j < ref_pos.size(); ++j) {
    //   if (i < target_pos.size() && ref_pos[j] == target_pos[i]) {
    //     lastValue = samples[i];
    //     ++i;
    //   }
    //   res[j] = lastValue;
    // }

    // return res;
    //
    std::vector<std::tuple<unsigned int, size_t, size_t>> intervals;
    size_t i = 0;
    size_t start_idx = 0;
    unsigned int lastValue = samples[0];

    for (size_t j = 0; j < ref_pos.size(); ++j) {
      if (i + 1 < target_pos.size() && ref_pos[j] >= target_pos[i + 1]) {
        intervals.emplace_back(lastValue, start_idx, j);
        start_idx = j + 1;
        lastValue = samples[++i];
      }
    }
    if (start_idx < ref_pos.size())
      intervals.emplace_back(lastValue, start_idx, ref_pos.size() - 1);

    return compactIntervals(intervals);
  }
  std::vector<std::tuple<unsigned int, size_t, size_t>> compactIntervals(
      const std::vector<std::tuple<unsigned int, size_t, size_t>> &intervals) {

    if (intervals.empty())
      return {};

    std::vector<std::tuple<unsigned int, size_t, size_t>> compacted;
    unsigned int current_value = std::get<0>(intervals[0]);
    size_t start_idx = std::get<1>(intervals[0]);
    size_t end_idx = std::get<2>(intervals[0]);

    for (size_t i = 1; i < intervals.size(); ++i) {
      unsigned int value = std::get<0>(intervals[i]);
      size_t new_start = std::get<1>(intervals[i]);
      size_t new_end = std::get<2>(intervals[i]);

      if (value == current_value && new_start == end_idx + 1) {
        end_idx = new_end;
      } else {
        compacted.emplace_back(current_value, start_idx, end_idx);
        current_value = value;
        start_idx = new_start;
        end_idx = new_end;
      }
    }

    compacted.emplace_back(current_value, start_idx, end_idx);
    return compacted;
  }

  std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
  cp(const std::vector<std::vector<unsigned int>> &matrix) {
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> result;

    for (size_t i = 0; i + 1 < matrix.size(); i += 2) {
      std::vector<std::pair<unsigned int, unsigned int>> pairs;

      size_t min_size = std::min(matrix[i].size(), matrix[i + 1].size());

      for (size_t j = 0; j < min_size; ++j) {
        pairs.emplace_back(matrix[i][j], matrix[i + 1][j]);
      }

      result.push_back(pairs);
    }

    return result;
  }
  std::vector<unsigned int> slicing(std::vector<unsigned int> const &v, int X,
                                    int Y) {
    std::vector<unsigned int> vector(v.begin() + X, v.begin() + Y + 1);

    return vector;
  }

  std::vector<unsigned int> gh(std::string f, std::string sample_name,
                               unsigned int haplotype,
                               std::unordered_set<unsigned int> ps, htsFile *fp,
                               bcf_hdr_t *hdr) {
    // htsFile *fp = bcf_open(f.c_str(), "r");

    // bcf_hdr_t *hdr = bcf_hdr_read(fp);
    int sample_idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, sample_name.c_str());
    bcf1_t *rec = bcf_init();
    int ngt, *gt_arr = NULL;
    int ret;
    std::vector<unsigned int> h;
    while ((ret = bcf_read(fp, hdr, rec)) >= 0) {
      // Unpack the genotypes
      bcf_unpack(rec, BCF_UN_ALL);

      int32_t *gt = NULL; // Pointer to genotype array
      int ngt = 0;        // Number of genotypes

      // Extract genotype data
      int ret = bcf_get_genotypes(hdr, rec, &gt, &ngt);
      if (ret < 0 || ngt <= sample_idx * 2 + 1) {
        free(gt);
        continue;
      }

      // Extract the desired haplotype
      int allele = gt[sample_idx * 2 + haplotype]; // 0-based index
      // std::cout << bcf_gt_allele(allele) << "\n";
      if (bcf_gt_is_missing(allele)) {
        h.push_back(-1); // Missing data
      } else {
        h.push_back(bcf_gt_allele(allele)); // Convert to allele index
      }

      // Free genotype memory
      free(gt);
    }

    // Cleanup
    free(gt_arr);
    bcf_destroy(rec);
    // bcf_hdr_destroy(hdr);
    // bcf_close(fp);
    return h;
  }

  /**
   * @brief function to compute queries with thresholds from a transposed tsv
   * or vcf file and output them on a file
   * @param filename queries file
   * @param out output file
   * @param extend_matches bool to extende mathc with rows values
   * @param verbose bool for extra prints
   * @param vcf bool to indicate if the file is a vcf
   */
  void query_match2(const char *filename, const char *out, bool verbose = false,
                    bool macs = false) {
    bool extend_matches = true;
    if (macs) {
      std::ifstream input_matrix(filename);
      std::ofstream out_match(out);
      if (input_matrix.is_open()) {
        std::string header1;
        std::string header2;
        std::string line;
        std::string garbage;
        std::string new_column;
        getline(input_matrix, line);
        getline(input_matrix, line);
        std::vector<std::string> queries_panel;
        while (getline(input_matrix, line) && !line.empty()) {
          std::istringstream is_col(line);
          is_col >> garbage;
          if (garbage == "TOTAL_SAMPLES:") {
            break;
          }
          is_col >> garbage >> garbage >> garbage >> new_column;
          queries_panel.push_back(new_column);
        }
        input_matrix.close();
        std::string query;
        std::vector<std::string> queries;
        if (out_match.is_open()) {
          for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
            if (verbose) {
              std::cout << i << ": \n";
            }
            for (auto &j : queries_panel) {
              query.push_back(j[i]);
            }
            queries.push_back(query);
            query.clear();
          }

          auto n_queries = queries.size();
          std::vector<ms_matches> matches_vec(n_queries);

#pragma omp parallel for default(none)                                         \
    shared(queries, matches_vec, n_queries, extend_matches, verbose)
          for (unsigned int i = 0; i < n_queries; i++) {
            matches_vec[i] =
                this->match_thr(queries[i], extend_matches, verbose);
          }
          if (extend_matches) {
            for (unsigned int i = 0; i < queries.size(); i++) {
              if (!matches_vec[i].haplos.empty()) {
                for (unsigned int j = 0;
                     j < matches_vec[i].basic_matches.size(); j++) {
                  auto len = std::get<1>(matches_vec[i].basic_matches[j]);
                  auto end = std::get<2>(matches_vec[i].basic_matches[j]);
                  for (unsigned int k = 0; k < matches_vec[i].haplos[j].size();
                       k++) {
                    out_match << "MATCH\t" << i << "\t"
                              << matches_vec[i].haplos[j][k] << "\t"
                              << end - (len - 1) << "\t" << end << "\t" << len
                              << "\n";
                  }
                }
              }
            }
          } else {
            for (unsigned int i = 0; i < queries.size(); i++) {
              for (unsigned int j = 0; j < matches_vec[i].basic_matches.size();
                   j++) {
                auto len = std::get<1>(matches_vec[i].basic_matches[j]);
                auto pos = std::get<0>(matches_vec[i].basic_matches[j]);
                auto end = std::get<2>(matches_vec[i].basic_matches[j]);
                out_match << "MATCH\t" << i << "\t" << pos << "\t"
                          << end - (len - 1) << "\t" << end << "\t" << len
                          << "\n";
              }
            }
          }
          out_match.close();
        } else {
          throw FileNotFoundException{};
        }

      } else {
        throw FileNotFoundException{};
      }
    } else {
      std::ofstream out_match(out);
      htsFile *fp = hts_open(filename, "rb");
      std::cout << "Reading VCF file...\n";

      bcf_hdr_t *hdr = bcf_hdr_read(fp);
      bcf1_t *rec = bcf_init();
      std::vector<std::string> queries_panel;
      std::string new_column;
      while (bcf_read(fp, hdr, rec) >= 0) {
        new_column = "";
        bcf_unpack(rec, BCF_UN_ALL);
        // read SAMPLE
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
        ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        int max_ploidy = ngt / nsmpl;
        for (i = 0; i < nsmpl; i++) {
          int32_t *ptr = gt_arr + i * max_ploidy;
          for (j = 0; j < max_ploidy; j++) {
            // if true, the sample has smaller ploidy
            if (ptr[j] == bcf_int32_vector_end)
              break;

            // missing allele
            if (bcf_gt_is_missing(ptr[j]))
              exit(-1);

            // the VCF 0-based allele index
            int allele_index = bcf_gt_allele(ptr[j]);
            new_column += std::to_string(allele_index);
          }
        }
        free(gt_arr);
        queries_panel.push_back(new_column);
      }
      bcf_hdr_destroy(hdr);
      hts_close(fp);
      bcf_destroy(rec);
      std::string query;
      std::vector<std::string> queries;
      if (out_match.is_open()) {
        for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
          if (verbose) {
            std::cout << i << ": \n";
          }
          for (auto &j : queries_panel) {
            query.push_back(j[i]);
          }
          queries.push_back(query);
          query.clear();
        }

        auto n_queries = queries.size();

        std::vector<ms_matches> matches_vec(n_queries);

#pragma omp parallel for default(none) shared(                                 \
        queries, matches_vec, n_queries, extend_matches, verbose, std::cout)
        for (unsigned int i = 0; i < n_queries; i++) {
          // std::cout <<  "query: " << i << "\n";
          matches_vec[i] = this->match_thr(queries[i], extend_matches, verbose);
        }

        if (extend_matches) {
          for (unsigned int i = 0; i < queries.size(); i++) {
            if (!matches_vec[i].haplos.empty()) {
              for (unsigned int j = 0; j < matches_vec[i].basic_matches.size();
                   j++) {
                auto len = std::get<1>(matches_vec[i].basic_matches[j]);
                auto end = std::get<2>(matches_vec[i].basic_matches[j]);
                for (unsigned int k = 0; k < matches_vec[i].haplos[j].size();
                     k++) {
                  out_match << "MATCH\t" << i << "\t"
                            << matches_vec[i].haplos[j][k] << "\t"
                            << end - (len - 1) << "\t" << end << "\t" << len
                            << "\n";
                }
              }
            }
          }
        } else {
          for (unsigned int i = 0; i < queries.size(); i++) {
            for (unsigned int j = 0; j < matches_vec[i].basic_matches.size();
                 j++) {
              auto len = std::get<1>(matches_vec[i].basic_matches[j]);
              auto pos = std::get<0>(matches_vec[i].basic_matches[j]);
              auto end = std::get<2>(matches_vec[i].basic_matches[j]);
              out_match << "MATCH\t" << i << "\t" << pos << "\t"
                        << end - (len - 1) << "\t" << end << "\t" << len
                        << "\n";
            }
          }
        }
        out_match.close();
      } else {
        throw FileNotFoundException{};
      }
    }
  }

  /**
   * function to get prefix array at a column
   * @param col required column
   * @return prefix array
   */
  std::vector<unsigned int> get_prefix(unsigned int col) {
    if (this->height == 0) {
      return {};
    }
    std::vector<unsigned int> pref;
    auto start_row = this->cols[col].sample_beg[0];
    pref.push_back(start_row);
    if (this->height == 1) {
      return pref;
    }
    auto next = this->phi->phi_inv(start_row, col);
    while (next.has_value()) {
      pref.push_back(next.value());
      next = this->phi->phi_inv(next.value(), col);
    }
    return pref;
  }

  /**
   * function to get divergence array at a column
   * @param col required column
   * @return divergence array
   */
  std::vector<unsigned int> get_divergence(unsigned int col) {
    if (this->height == 0) {
      return {};
    }
    std::vector<unsigned int> div;
    auto start_row = this->cols[col].sample_beg[0];
    div.push_back(this->phi->plcp(start_row, col));
    if (this->height == 1) {
      return div;
    }
    auto next = this->phi->phi_inv(start_row, col);
    while (next.has_value()) {
      div.push_back(this->phi->plcp(next.value(), col));
      next = this->phi->phi_inv(next.value(), col);
    }
    return div;
  }

  /**
   * function to get prefix/divergence array at a column
   * @param col required column
   * @return prefix and divergence array
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  get_prefix_divergence(unsigned int col) {
    if (this->height == 0) {
      return {};
    }
    std::vector<std::pair<unsigned int, unsigned int>> prefdiv;
    auto start_row = this->cols[col].sample_beg[0];
    prefdiv.push_back(
        std::make_pair(start_row, this->phi->plcp(start_row, col)));
    if (this->height == 1) {
      return prefdiv;
    }
    auto next = this->phi->phi_inv(start_row, col);
    while (next.has_value()) {
      prefdiv.push_back(
          std::make_pair(next.value(), this->phi->plcp(next.value(), col)));
      next = this->phi->phi_inv(next.value(), col);
    }
    return prefdiv;
  }

  /**
   * function to get u/v arrays at a column
   * @param col required column
   * @return u/v arrays
   */
  std::vector<std::pair<unsigned int, unsigned int>> get_u_v(unsigned int col) {
    if (this->height == 0) {
      return {};
    }
    std::vector<std::pair<unsigned int, unsigned int>> uv;
    for (unsigned int i = 0; i < this->height; i++) {
      auto tmp = uvtrick(col, index_to_run(i, col));
      uv.push_back(tmp);
    }
    return uv;
  }

  /**
   * function to get a column of the pbwt
   * @param col required column
   * @return column as string
   */
  std::string get_col(unsigned int col) {
    if (this->height == 1) {
      return {};
    }
    std::string c = "";
    bool s = this->cols[col].zero_first;
    for (unsigned int i = 0; i < this->height; i++) {
      c += get_next_char(s, index_to_run(i, col));
    }
    return c;
  }

  std::string get_orig_col(unsigned int col) {
    auto pbwt_col = get_col(col);
    auto p = get_prefix(col);
    std::string orig_col(pbwt_col.size(), '0');
    for (unsigned int i = 0; i < pbwt_col.size(); i++) {
      if (pbwt_col[i] == '1')
        orig_col[p[i]] = '1';
    }
    return orig_col;
  }

  std::vector<std::vector<unsigned int>>
  get_sub_matrix(std::vector<unsigned int> rows, unsigned int cb, unsigned ce) {
    std::vector<std::vector<unsigned int>> sub_m(
        rows.size(), std::vector<unsigned int>(ce - cb + 1));
    std::vector<std::string> cols;
    for (unsigned int i = cb; i <= ce; i++)
      cols.push_back(get_orig_col(i));
    unsigned int i = 0;
    for (auto col : cols) {
      for (unsigned int j = 0; j < rows.size(); j++) {
        sub_m[j][i] = (unsigned int)(col[rows[j]] - '0');
      }
      i++;
    }
    return sub_m;
  }

  bool check_comp(const std::vector<unsigned int> &row1,
                  const std::vector<unsigned int> &row2) {
    if (row1.size() != row2.size())
      return false;
    for (size_t i = 0; i < row1.size(); ++i) {
      if (row1[i] == row2[i])
        return false;
    }
    return true;
  }

  std::vector<std::pair<unsigned int, unsigned int>>
  get_comp(const std::vector<std::vector<unsigned int>> &mat,
           std::vector<unsigned int> r) {
    std::vector<std::pair<unsigned int, unsigned int>> result;
    // #pragma omp parallel for
    for (size_t i = 0; i < mat.size(); ++i) {
      // #pragma omp parallel for
      for (size_t j = i + 1; j < mat.size(); ++j) {
        if (check_comp(mat[i], mat[j])) {
          if (r[i] < r[j])
            result.emplace_back(r[i], r[j]);
          else
            result.emplace_back(r[j], r[i]);
        }
      }
    }
    return result;
  }
  /**
   * function to get the total number of runs in the RLPBWT
   * @return total number of run
   */
  unsigned int get_run_number() {
    unsigned int count_run = 0;
    for (unsigned int i = 0; i < this->cols.size(); ++i) {
      count_run += cols[i].sample_beg.size();
    }
    return count_run;
  }

  /**
   * function to print in runs.txt the number of run in every
   * column
   */
  void get_run_col(char *filename) {
    std::ofstream myfile;
    myfile.open(filename);
    unsigned int min = 1;
    unsigned int max = 1;
    for (unsigned int i = 0; i < this->cols.size(); ++i) {
      auto x = cols[i].sample_beg.size();
      myfile << x << " ";
      if (x < min) {
        min = x;
      }
      if (x > max) {
        max = x;
      }
    }
    myfile << "\n" << min << " " << max;
  }

  /**
   * function to get the total number of phi/phi_inv element() in the RLPBWT
   * @return total number of run
   */
  std::pair<unsigned int, unsigned int> get_phi_number() {
    if (this->phi) {
      unsigned int count_phi = 0;
      for (unsigned int i = 0; i < this->phi->phi_supp.size(); ++i) {
        count_phi += this->phi->phi_supp[i].size();
      }
      unsigned int count_phi_inv = 0;
      for (unsigned int i = 0; i < this->phi->phi_inv_supp.size(); ++i) {
        count_phi_inv += this->phi->phi_inv_supp[i].size();
      }
      return std::make_pair(count_phi, count_phi_inv);
    }
  }

  /**
   * @brief function to obtain size in bytes of the matching statistics
   * supported RLPBWT
   * @param verbose bool for extra prints
   * @return size in bytes
   */
  unsigned long long size_in_bytes(bool verbose = false) {
    unsigned long long size = 0;
    unsigned long long size_run = 0;
    unsigned long long size_thr = 0;
    unsigned long long size_uv = 0;
    auto lp_size = sdsl::size_in_bytes(this->last_pref);
    unsigned long long size_samples = lp_size;
    size += lp_size;
    for (unsigned int i = 0; i < this->cols.size(); ++i) {
      size += this->cols[i].size_in_bytes();
      size_run += sdsl::size_in_bytes(this->cols[i].p);
      size_thr += sdsl::size_in_bytes(this->cols[i].t);
      size_uv += sdsl::size_in_bytes(this->cols[i].uv);
      size_samples += sdsl::size_in_bytes(this->cols[i].sample_beg) +
                      sdsl::size_in_bytes(this->cols[i].sample_end);
    }
    if (verbose) {
      std::cout << "run: " << size_run << " bytes\n";
      std::cout << "thr: " << size_thr << " bytes\n";

      std::cout << "uv: " << size_uv << " bytes\n";
      std::cout << "samples: " << size_samples << " bytes\n";
      std::cout << "rlpbwt (with also c values and other support variables): "
                << size << " bytes\n";
    }
    size += (sizeof(bool) * 2);
    size += (sizeof(unsigned int) * 2);
    if (this->is_extended) {
      auto size_phi = this->phi->size_in_bytes(verbose);
      size += size_phi;
      // std::cout << "phi support: " << size_phi << " bytes\n";
    }
    return size;
  }

  /**
   * @brief function to obtain size in megabytes of the matching statistics
   * supported RLPBWT
   * @param verbose bool for extra prints
   * @return size in megabytes
   */
  double size_in_mega_bytes(bool verbose = true) {
    double size = 0;
    double to_mega = ((double)1 / (double)1024) / (double)1024;
    double size_run = 0;
    double size_thr = 0;
    double size_uv = 0;
    auto lp_size = sdsl::size_in_mega_bytes(this->last_pref);
    double size_samples = lp_size;
    size += lp_size;
    for (unsigned int i = 0; i < this->cols.size(); ++i) {
      size += this->cols[i].size_in_mega_bytes();
      size_run += sdsl::size_in_mega_bytes(this->cols[i].p);
      size_thr += sdsl::size_in_mega_bytes(this->cols[i].t);
      size_uv += sdsl::size_in_mega_bytes(this->cols[i].uv);
      size_samples += sdsl::size_in_mega_bytes(this->cols[i].sample_beg) +
                      sdsl::size_in_mega_bytes(this->cols[i].sample_end) +
                      sdsl::size_in_mega_bytes(this->cols[i].sample_beg_lcp);
    }
    size += (sizeof(bool) * 2 * to_mega);
    size += (sizeof(unsigned int) * 2 * to_mega);
    if (verbose) {
      std::cout << "run: " << size_run << " megabytes\n";
      std::cout << "thr: " << size_thr << " megabytes\n";
      std::cout << "uv: " << size_uv << " megabytes\n";
      std::cout << "samples: " << size_samples << " megabytes\n";
      std::cout << "rlpbwt (mapping): " << size << " megabytes\n";
    }

    if (this->is_extended) {
      auto size_phi = this->phi->size_in_mega_bytes(verbose);
      size += size_phi;
    }
    return size;
  }

  /**
   * @brief function to serialize the matching statistics supported RLPBWT
   * @param out std::ostream object to stream the serialization
   * @return size of the serialization
   */
  size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                   const std::string &name = "") {
    sdsl::structure_tree_node *child =
        sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;
    out.write((char *)&this->width, sizeof(this->width));
    written_bytes += sizeof(this->width);
    out.write((char *)&this->height, sizeof(this->height));
    written_bytes += sizeof(this->height);
    out.write((char *)&this->k_smem, sizeof(this->k_smem));
    written_bytes += sizeof(this->k_smem);
    out.write((char *)&this->is_extended, sizeof(this->is_extended));
    written_bytes += sizeof(this->is_extended);
    for (unsigned int i = 0; i < this->cols.size(); i++) {
      std::string label = "col_" + std::to_string(i);
      written_bytes += this->cols[i].serialize(out, child, label);
    }
    written_bytes += this->last_pref.serialize(out, child, "last_pref");
    written_bytes += this->last_div.serialize(out, child, "last_div");
    if (this->is_extended) {
      written_bytes += this->phi->serialize(out, child, "phi");
    }
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  /**
   * @brief function to load the matching statistics supported RLPBWT
   * object
   * @param in std::istream object from which load the matching statistics
   * supported RLPBWT structure object
   */
  void load(std::istream &in) {
    in.read((char *)&this->width, sizeof(this->width));
    in.read((char *)&this->height, sizeof(this->height));
    in.read((char *)&this->k_smem, sizeof(this->k_smem));
    in.read((char *)&this->is_extended, sizeof(this->is_extended));
    this->cols = std::vector<rl_column>(this->width + 1);
    for (unsigned int i = 0; i <= this->width; i++) {
      this->cols[i].load(in);
    }
    this->last_pref.load(in);
    this->last_div.load(in);
    if (this->is_extended) {
      this->phi = new phi_ds();
      this->phi->load(in);
    }
  }
};

#endif // RLPBWT_RLPBWT_NAIVE_MS_H
