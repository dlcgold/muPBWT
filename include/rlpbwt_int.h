//
// Created by dlcgold on 16/08/22.
//

#ifndef RLPBWT_RLPBWT_NAIVE_MS_H
#define RLPBWT_RLPBWT_NAIVE_MS_H

#include <type_traits>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "rl_column.h"
#include "exceptions.h"
#include "htslib/vcf.h"
#include "ms.h"
#include "ms_matches.h"
#include "phi_ds.h"


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
    std::vector <rl_column> cols;

    /**
     * @brief phi/phi_inv support data structure
     */
    phi_ds *phi = nullptr;

private:
    /**
     * @brief compressed int vector for last prefix array
     */
    sdsl::int_vector<> last_pref;
    sdsl::int_vector<> last_div;

    /**
     * @brief function to build a naive column for a column of the panel
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     * @return the new naive column
     */
//    rl_column build_column(std::string &column,
//                           std::vector<unsigned int> &pref,
//                           sdsl::int_vector<> &div) {
    rl_column build_column(std::string &column,
                           std::vector<unsigned int> &pref,
                           std::vector<unsigned int> &div,
                           std::vector <intv> &supp_b,
                           std::vector <intv> &supp_e,
                           unsigned int col_i = 0) {
        //unsigned int height = pref.size();
        // variable for "c" value
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
        std::vector <std::pair<unsigned int, unsigned int>> rows;
        std::vector<unsigned int> thr;
        // support vector to store prefix array samples
        std::vector <std::pair<unsigned int, unsigned int>> samples;
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
            if ((i == height - 1) ||
                (column[pref[i]] != column[pref[i + 1]])) {
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

        sdsl::int_vector<>
                t_vec(rows.size());
        sdsl::int_vector<>
                p_vec(rows.size());
        sdsl::int_vector<>
                uv_vec(rows.size());
        sdsl::int_vector<>
                sb_vec(rows.size());
        sdsl::int_vector<>
                se_vec(rows.size());
        sdsl::int_vector<>
                sbl_vec(rows.size());
        sdsl::int_vector<>
                sel_vec(rows.size());
        for (unsigned int i = 0; i < rows.size(); i++) {
            p_vec[i] = rows[i].first;
            uv_vec[i] = rows[i].second;
            t_vec[i] = thr[i];
            sb_vec[i] = samples[i].first;
            se_vec[i] = samples[i].second;
            sbl_vec[i] = samples_lcp[i];
            supp_b[samples[i].first].push_back(col_i);
            supp_e[samples[i].second].push_back(col_i);
            //supp_l[samples_lcp[i]].push_back(count);
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
        //sdsl::int_vector<> new_div(height);
        std::vector<unsigned int> new_div(height);
        unsigned int count0 = 0;
        unsigned int lcs = -1;

        for (unsigned int i = 0; i < height; i++) {
            //lcs = std::min(lcs, static_cast<unsigned int>(div[i]));
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
            //lcs = std::min(lcs, static_cast<unsigned int>(div[i]));
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
             get_next_char(this->cols[col_index].zero_first,
                           run_index) == '1') ||
            (symbol == '1' &&
             get_next_char(this->cols[col_index].zero_first,
                           run_index) == '0')) {
            offset = 0;
        }
        // obtain "u" and "v"
        auto uv = uvtrick(col_index, run_index);
        if (verbose) {
            std::cout << uv.first << ", " << uv.second << "\n";
        }
        // fix for the last index that's "outside" the column
        if (this->cols[col_index].p[run_index] + offset == this->height) {
            if (get_next_char(this->cols[col_index].zero_first,
                              run_index) ==
                '0') {
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

    unsigned int
    reverse_lf(unsigned int col_index, unsigned int index, bool verbose) const {
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
        // - if index is less than previous "c" it means that it comes from a zero
        //   element, so we will search the correct u value
        // - otherwise it means that it comes from a one element, so we will search
        //   the correct v value
        if (index < c) {
            u = index;
            if (verbose) {
                std::cout << "u: " << u << "\n";
            }
            // iteration over the u values to find the correct one
            unsigned int prevu = 0;
            unsigned int nextu = 0;
            for (unsigned int i = 0;
                 i < this->cols[col_index].p.size() - 1; i++) {
                prevu = uvtrick(col_index, i).first;
                nextu = uvtrick(col_index, i + 1).first;
                if (prevu <= u && u < nextu) {
                    pos = i;
                    found = true;
                    break;
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
            unsigned int prevv = 0;
            unsigned int nextv = 0;
            v = index - c;
            if (verbose) {
                std::cout << "v: " << v << "\n";
            }
            // iteration over the v values to find the correct one
            for (unsigned int i = 0;
                 i < this->cols[col_index].p.size() - 1; ++i) {
                prevv = uvtrick(col_index, i).second;
                nextv = uvtrick(col_index, i + 1).second;
                if (prevv <= v && v < nextv) {
                    pos = i;
                    found = true;
                    break;
                }
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
    unsigned int index_to_run(unsigned int index,
                              unsigned int col_index) const {
        // if requested index is equal or greater than the p value of the last
        // run return the index of the last run
        if (index >=
            this->cols[col_index].p[this->cols[col_index].p.size() - 1]) {
            return this->cols[col_index].p.size() - 1;
        }

        // binary search to compute run index
        unsigned int bi = 0;
        unsigned int e = this->cols[col_index].p.size();
        unsigned int pos = (e - bi) / 2;
        while (pos != e && this->cols[col_index].p[pos] != index) {
            if (index < (unsigned int) this->cols[col_index].p[pos]) {
                e = pos;
            } else {
                if (pos + 1 == e ||
                    (unsigned int) this->cols[col_index].p[pos + 1] >
                    index) {
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
    std::pair<unsigned int, unsigned int> uvtrick(
            unsigned int col_index, unsigned int run_index) const {
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
    void
    extend_haplos(ms_matches &ms_matches, std::vector<unsigned int> ms_supp) {

        // iterate over every basic match
        for (unsigned int i = 0;
             i < ms_matches.basic_matches.size(); i++) {
            // initialize the vector that will contain row indices
            std::vector<unsigned int> haplos;

            // extract information from the current basic match
            unsigned int start_row = std::get<0>(
                    ms_matches.basic_matches[i]);
            haplos.emplace_back(start_row);
            unsigned int curr_len = std::get<1>(
                    ms_matches.basic_matches[i]);
            unsigned int curr_col = std::get<2>(
                    ms_matches.basic_matches[i]);
            //std::cout << "for match " << start_row << " " << curr_len << "\n";
            // initialize boolean and temporary variables for go up/down in
            // search of matching rows
            bool check_down = true;
            unsigned int down_row = 0;
            bool check_up = true;
            unsigned int up_row = 0;

            char curr_s = get_next_char(this->cols[curr_col].zero_first,
                                        index_to_run(ms_supp[curr_col],
                                                     curr_col));
            unsigned int down_index =
                    lf(curr_col, ms_supp[curr_col], curr_s) + 1;
            unsigned int up_index = lf(curr_col, ms_supp[curr_col], curr_s);
            // go down/up and add row that has a lce of at least current len in
            // matching statistics (if panel is not saved as SLP the computation
            // of lce is simulated with random access to the panel)

            // down
            while (check_down) {
                auto phi_res = this->phi->phi_inv(start_row,
                                                  curr_col + 1);

                if (!phi_res.has_value()) {
                    break;
                }
                down_row = phi_res.value();
                //std::cout << "check d: " << down_row << " " << down_index << " " << curr_col << " " << start_row << "\n";
                //std::cout << curr_len << " "
                //          << this->phi->plcp(down_row, curr_col) << "\n";
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
                auto phi_res = this->phi->phi(start_row,
                                              curr_col + 1);

                if (!phi_res.has_value()) {
                    break;
                }
                up_row = phi_res.value();
                //std::cout << "check u: " << up_row << " " << up_index << " " << curr_col << " " << start_row << "\n";
                //std::cout << curr_len << " "
                //          << this->phi->plcp(start_row, curr_col) << "\n";
//                if (up_index < i_r || up_index > e_r) {
//                    break;
//                }

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
            //ref.clear();
        }
    }

//    bool check_extend(unsigned int col, std::vector<char> f_r,
//                      unsigned int s_r,
//                      unsigned int length) const {
//        unsigned int start = 0;
//        if ((int) col - ((int) length - 1) >= 0) {
//            start = col - (length - 1);
//        }
//        auto size = f_r.size();
//        bool check = true;
//        unsigned int end = f_r.size() - 1;
//        s_r = reverse_lf(col + 1, s_r, false);
//        for (unsigned int i = col; i >= start; i--) {
//            if (i >= 0) {
//                char s = get_next_char(this->cols[i].zero_first,
//                                       index_to_run(s_r, i));
//                if (f_r[end] != s) {
//                    check = false;
//                    break;
//                }
//
//                s_r = reverse_lf(i, s_r, false);
//
//                end--;
//            }
//            if (i == 0) {
//                break;
//            }
//        }
//        return check;
//    }


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
    explicit rlpbwt_int(const char *filename,
                        bool verbose = false,
                        bool macs = false) {
        if (macs) {
            std::ifstream input_matrix(filename);
            if (input_matrix.is_open()) {
                std::string header1;
                std::string header2;
                std::string line;
                std::string garbage;
                std::string new_column;

                getline(input_matrix, header1);
                getline(input_matrix, header2);
                getline(input_matrix, line);

                std::istringstream is(line);
                is >> garbage >> garbage >> garbage >> garbage
                   >> new_column;
                unsigned int tmp_height = new_column.size();
                std::cout << "Total samples: " << tmp_height << "\n";
                /*unsigned int tmp_width = std::count(
                        std::istreambuf_iterator<char>(input_matrix),
                        std::istreambuf_iterator<char>(), '\n');
                        */
                //auto tmp_width = 1;
//                while (getline(input_matrix, line) && !line.empty()) {
//                    std::istringstream is_col(line);
//                    is_col >> garbage;
//                    if (garbage == "TOTAL_SAMPLES:") {
//                        break;
//                    }
//                    tmp_width++;
//                }
//                std::cout << "w: " << tmp_width << "\n";
                this->width = 0;
                this->height = tmp_height;
                input_matrix.clear();
                input_matrix.seekg(0, std::ios::beg);
                //this->cols = std::vector<rl_column>(tmp_width + 1);
                std::vector<unsigned int> pref(tmp_height);
                //sdsl::int_vector<> div(tmp_height);
                std::vector<unsigned int> div(tmp_height);
                this->last_pref.resize(tmp_height);
                this->last_div.resize(tmp_height);
                auto supp_b = std::vector<intv>(this->height);
                auto supp_e = std::vector<intv>(this->height);
                for (unsigned int i = 0; i < tmp_height; i++) {
                    pref[i] = i;
                    div[i] = 0;
                }
                unsigned int count = 0;
                std::string last_col;
                getline(input_matrix, line);
                getline(input_matrix, line);
                std::string last_column;
                while (getline(input_matrix, line) && !line.empty()) {
                    std::cout << count << "\r";
                    std::istringstream is_col(line);
                    is_col >> garbage;
                    if (garbage == "TOTAL_SAMPLES:") {
                        break;
                    }
                    is_col >> garbage >> garbage >> garbage >> new_column;
                    if (verbose) {
                        std::cout << "\nnew_column " << count << "\n";
                        std::cout << new_column << "\n"
                                  << this->cols[count]
                                  << "\n-------------------------------\n";
                    }
                    auto col = rlpbwt_int::build_column(new_column, pref,
                                                        div, supp_b, supp_e,
                                                        count);
                    //this->cols[count] = col;
                    this->cols.emplace_back(col);
                    rlpbwt_int::update(new_column, pref, div);
                    last_col = new_column;
                    count++;
                }
                for (unsigned int i = 0; i < pref.size(); i++) {
                    this->last_pref[i] = pref[i];
                    this->last_div[i] = div[i];
                }
                std::cout << std::endl;
                this->width = cols.size();
                std::cout << "Total sites: " << this->width << "\n";
                auto col = rlpbwt_int::build_column(last_col, pref, div, supp_b,
                                                    supp_e, count);
                //this->cols[count] = col;
                this->cols.emplace_back(col);
                for (unsigned int i = 0; i < this->height; i++) {
                    if (supp_b[i].v.size() == 0 ||
                        supp_b[i].v[supp_b[i].v.size() - 1] != this->width) {
                        supp_b[i].push_back(this->width);
                    }
                    supp_b[i].compress();
                    if (supp_e[i].v.size() == 0 ||
                        supp_e[i].v[supp_e[i].v.size() - 1] != this->width) {
                        supp_e[i].push_back(this->width);
                    }
                    supp_e[i].compress();
                }
                sdsl::util::bit_compress(this->last_pref);
                sdsl::util::bit_compress(this->last_div);
                this->phi = new
                        phi_ds(this->cols, this->height, this->width,
                               this->last_pref, this->last_div, supp_b, supp_e,
                               verbose);
                this->is_extended = true;
                input_matrix.close();
            } else {
                throw FileNotFoundException{};
            }
        } else {
            htsFile *fp = hts_open(filename, "rb");
            std::cout << "Reading VCF file...\n";
            if (fp == NULL) {
                throw FileNotFoundException{};
            }

            bcf_hdr_t *hdr = bcf_hdr_read(fp);
            bcf1_t *rec = bcf_init();
            this->height = bcf_hdr_nsamples(hdr) * 2;
            this->width = 0;
            //while (bcf_read(fp, hdr, rec) >= 0) {
            //    this->width++;
            //}
            std::cout << "Total haplotypes: " << this->height << "\n";
            //std::cout << "h: " << this->height << "\n";
            //std::cout << "w: " << this->width << "\n";
            // bcf_hdr_destroy(hdr);
            //hts_close(fp);
            // bcf_destroy(rec);
            fp = hts_open(filename, "rb");
            hdr = bcf_hdr_read(fp);
            rec = bcf_init();

            //this->cols = std::vector<rl_column>(this->width + 1);
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
                this->width++;
                std::cout << "Analyze site: " << count << "\r";
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
                        if (ptr[j] == bcf_int32_vector_end) break;

                        // missing allele
                        if (bcf_gt_is_missing(ptr[j])) exit(-1);

                        // the VCF 0-based allele index
                        int allele_index = bcf_gt_allele(ptr[j]);
                        new_column += std::to_string(allele_index);
                    }
                }
                free(gt_arr);
                auto col = rlpbwt_int::build_column(new_column, pref,
                                                    div, supp_b, supp_e, count);
                //this->cols[count] = col;
                this->cols.emplace_back(col);
                rlpbwt_int::update(new_column, pref, div);
                last_col = new_column;

                count++;
            }
            std::cout << std::endl;
            std::cout << "Total sites: " << this->width << "\n";
            for (unsigned int i = 0; i < pref.size(); i++) {
                this->last_pref[i] = pref[i];
                this->last_div[i] = div[i];
            }

            auto col = rlpbwt_int::build_column(last_col, pref, div, supp_b,
                                                supp_e, count);
            this->cols.emplace_back(col);
            for (unsigned int i = 0; i < this->height; i++) {
                if (supp_b[i].v.size() == 0 ||
                    supp_b[i].v[supp_b[i].v.size() - 1] != this->width) {
                    supp_b[i].push_back(this->width);
                }
                supp_b[i].compress();
                if (supp_e[i].v.size() == 0 ||
                    supp_e[i].v[supp_e[i].v.size() - 1] != this->width) {
                    supp_e[i].push_back(this->width);
                }
                supp_e[i].compress();
            }
            sdsl::util::bit_compress(this->last_pref);
            sdsl::util::bit_compress(this->last_div);
            this->phi = new phi_ds(this->cols, this->height, this->width,
                                   this->last_pref, this->last_div, supp_b,
                                   supp_e,
                                   verbose);
            this->is_extended = true;

            bcf_hdr_destroy(hdr);
            hts_close(fp);
            bcf_destroy(rec);
        }
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


    /**
     * @brief function to compute matching statistics matches with a given query
     * using thresholds
     * @param query haplotype query as std::string
     * @param extend_matches bool to check if extend matching statistics matches
     * with rows
     * @param verbose bool for extra prints
     * @return matching statistics matches
     * @attention use this function is enabled iff thresholds are calculated
     */
    ms_matches match_thr(const std::string &query, bool extend_matches = false,
                         bool verbose = false) {

        // compute the match iff |query| is equal to the width of the panel
        if (query.size() != this->width) {
            std::cout << query.size() << " != " << this->width
                      << "\n";
            throw NotEqualLengthException{};
        }
        // if required extend with the phi support struct (iff not already
        // extended)
        if (extend_matches && !this->is_extended) {
            //this->extend();
        }
        // initialize matching statistics
        ms ms(query.size());
        std::vector<unsigned int> ms_supp(query.size());

        // algorithm begin from the last row of the first column
        // so we obtain the prefix array value (from the samples), the run index
        // and the relative symbol
        auto curr_pos = static_cast<unsigned int>(
                this->cols[0].sample_end[this->cols[0].sample_end.size() -
                                         1]);
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
                std::cout << curr_index << " " << curr_run << " "
                          << curr_pos
                          << " " << symbol << "\n";
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
                    symbol =
                            get_next_char(this->cols[i + 1].zero_first,
                                          curr_run);
                    if (verbose) {
                        std::cout << "new: " << curr_index << " "
                                  << curr_run
                                  << " " << curr_pos << " " << symbol
                                  << "\n";
                    }
                }
            } else {
                // get threshold
                auto thr = this->cols[i].t[curr_run];

                if (this->cols[i].sample_beg.size() == 1) {
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
                                this->cols[i + 1].sample_end
                                [this->cols[i + 1].sample_end.size() - 1]);
                        curr_index = this->height - 1;

                        curr_run = index_to_run(curr_index, i + 1);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos << " " << symbol << "\n";
                        }
                    }
                } else if (curr_run != 0 &&
                           ((curr_index < thr) ||
                            curr_run ==
                            this->cols[i].sample_beg.size() - 1)) {
                    // if we are above the threshold we go up (if we are not in
                    // the first run). We also go up if we are in the last run
                    if (verbose) {
                        std::cout << "mismatch_up: ";
                    }
                    curr_index = this->cols[i].p[curr_run] - 1;
                    curr_pos = static_cast<unsigned int>(
                            this->cols[i].sample_end[curr_run - 1]);
                    if (verbose) {
                        std::cout << "update: " << curr_index << " "
                                  << curr_pos
                                  << " " << symbol << "\n";
                    }
                    // report in matching statistics row vector
                    ms.row[i] = curr_pos;

                    ms_supp[i] = curr_index;

                    // update index, run, symbol if we are not at the end
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = index_to_run(curr_index, i + 1);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "new: " << curr_index << " "
                                      << curr_run << " " << curr_pos << " "
                                      << symbol << "\n";
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
                    ms_supp[i] = curr_index;

                    if (verbose) {
                        std::cout << "update: " << curr_index << " "
                                  << curr_pos
                                  << " " << symbol << "\n";
                    }
                    // update index, run, symbol if we are not at the end
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = index_to_run(curr_index, i + 1);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "new: " << curr_index << " "
                                      << curr_run << " " << curr_pos << " "
                                      << symbol << "\n";
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
            } else if (i != 0 && ms.row[i] == ms.row[i - 1] &&
                       ms.len[i - 1] != 0) {
                ms.len[i] = ms.len[i - 1] + 1;
            } else {

                int tmp_index = (int) i;
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
        // initialize struct for matches
        ms_matches ms_matches;
        // save every match from matching statistics (when we have a "peak" in
        // ms len vector)
        for (unsigned int i = 0; i < ms.len.size(); i++) {
            if ((i != ms.len.size() - 1 && ms.len[i] > 0 &&
                 ms.len[i] >= ms.len[i + 1]) ||
                (i == ms.len.size() - 1 && ms.len[i] != 0)) {
                ms_matches.basic_matches.emplace_back(ms.row[i], ms.len[i],
                                                      i);
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
    }


    /**
     * @brief function to compute queries with thresholds from a transposed tsv or vcf
     * file and output them on a file
     * @param filename queries file
     * @param out output file
     * @param extend_matches bool to extende mathc with rows values
     * @param verbose bool for extra prints
     * @param vcf bool to indicate if the file is a vcf
     */
    void query_match(const char *filename, const char *out,
                     bool verbose = false,
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
                std::vector <std::string> queries_panel;
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
                std::vector <std::string> queries;
                if (out_match.is_open()) {
                    for (unsigned int i = 0;
                         i < queries_panel[0].size(); i++) {
                        if (verbose) {
                            std::cout << i << ": \n";
                        }
                        for (auto &j: queries_panel) {
                            query.push_back(j[i]);
                        }
                        queries.push_back(query);
                        query.clear();
                    }

                    auto n_queries = queries.size();
                    std::vector <ms_matches> matches_vec(n_queries);
#pragma omp parallel for default(none) \
    shared(queries, matches_vec, n_queries, extend_matches, verbose)
                    for (unsigned int i = 0; i < n_queries; i++) {
                        // std::cout << i << "\n";
                        matches_vec[i] = this->match_thr(
                                queries[i], extend_matches, verbose);
                    }
                    if (extend_matches) {
                        for (unsigned int i = 0; i < queries.size(); i++) {
                            if (!matches_vec[i].haplos.empty()) {
                                for (unsigned int j = 0;
                                     j <
                                     matches_vec[i].basic_matches.size();
                                     j++) {
                                    auto len = std::get<1>(
                                            matches_vec[i].basic_matches[j]);
                                    auto end = std::get<2>(
                                            matches_vec[i].basic_matches[j]);
                                    for (unsigned int k = 0;
                                         k <
                                         matches_vec[i].haplos[j].size();
                                         k++) {
                                        out_match << "MATCH\t" << i << "\t"
                                                  << matches_vec[i].haplos[j][k]
                                                  << "\t" << end - (len - 1)
                                                  << "\t" << end << "\t"
                                                  << len
                                                  << "\n";
                                    }
                                }
                            }
                        }
                    } else {
                        for (unsigned int i = 0; i < queries.size(); i++) {
                            for (unsigned int j = 0;
                                 j <
                                 matches_vec[i].basic_matches.size(); j++) {
                                auto len = std::get<1>(
                                        matches_vec[i].basic_matches[j]);
                                auto pos = std::get<0>(
                                        matches_vec[i].basic_matches[j]);
                                auto end = std::get<2>(
                                        matches_vec[i].basic_matches[j]);
                                out_match << "MATCH\t" << i << "\t" << pos
                                          << "\t" << end - (len - 1) << "\t"
                                          << end << "\t" << len << "\n";
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
            std::vector <std::string> queries_panel;
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
                        if (ptr[j] == bcf_int32_vector_end) break;

                        // missing allele
                        if (bcf_gt_is_missing(ptr[j])) exit(-1);

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
            std::vector <std::string> queries;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    if (verbose) {
                        std::cout << i << ": \n";
                    }
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    queries.push_back(query);
                    query.clear();
                }

                auto n_queries = queries.size();

                std::vector <ms_matches> matches_vec(n_queries);

#pragma omp parallel for default(none) \
    shared(queries, matches_vec, n_queries, extend_matches, verbose)
                for (unsigned int i = 0; i < n_queries; i++) {
                    matches_vec[i] =
                            this->match_thr(queries[i], extend_matches,
                                            verbose);
                }

                if (extend_matches) {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        if (!matches_vec[i].haplos.empty()) {
                            for (unsigned int j = 0;
                                 j <
                                 matches_vec[i].basic_matches.size(); j++) {
                                auto len = std::get<1>(
                                        matches_vec[i].basic_matches[j]);
                                auto end = std::get<2>(
                                        matches_vec[i].basic_matches[j]);
                                for (unsigned int k = 0;
                                     k <
                                     matches_vec[i].haplos[j].size(); k++) {
                                    out_match << "MATCH\t" << i << "\t"
                                              << matches_vec[i].haplos[j][k]
                                              << "\t" << end - (len - 1)
                                              << "\t"
                                              << end << "\t" << len << "\n";
                                }
                            }
                        }
                    }
                } else {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        for (unsigned int j = 0;
                             j < matches_vec[i].basic_matches.size(); j++) {
                            auto len =
                                    std::get<1>(
                                            matches_vec[i].basic_matches[j]);
                            auto pos =
                                    std::get<0>(
                                            matches_vec[i].basic_matches[j]);
                            auto end =
                                    std::get<2>(
                                            matches_vec[i].basic_matches[j]);
                            out_match << "MATCH\t" << i << "\t" << pos
                                      << "\t"
                                      << end - (len - 1) << "\t" << end
                                      << "\t"
                                      << len << "\n";
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
    std::vector <std::pair<unsigned int, unsigned int>>
    get_prefix_divergence(unsigned int col) {
        if (this->height == 0) {
            return {};
        }
        std::vector <std::pair<unsigned int, unsigned int>> prefdiv;
        auto start_row = this->cols[col].sample_beg[0];
        prefdiv.push_back(
                std::make_pair(start_row, this->phi->plcp(start_row, col)));
        if (this->height == 1) {
            return prefdiv;
        }
        auto next = this->phi->phi_inv(start_row, col);
        while (next.has_value()) {
            prefdiv.push_back(std::make_pair(next.value(),
                                             this->phi->plcp(next.value(),
                                                             col)));
            next = this->phi->phi_inv(next.value(), col);
        }
        return prefdiv;
    }

    /**
   * function to get u/v arrays at a column
   * @param col required column
   * @return u/v arrays
   */
    std::vector <std::pair<unsigned int, unsigned int>>
    get_u_v(unsigned int col) {
        if (this->height == 0) {
            return {};
        }
        std::vector <std::pair<unsigned int, unsigned int>> uv;
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
            for (unsigned int i = 0;
                 i < this->phi->phi_inv_supp.size(); ++i) {
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
            std::cout
                    << "rlpbwt (with also c values and other support variables): "
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
        double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
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
            size_samples +=
                    sdsl::size_in_mega_bytes(this->cols[i].sample_beg) +
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
            std::cout
                    << "rlpbwt (mapping): "
                    << size << " megabytes\n";
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
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(
                v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        out.write((char *) &this->width, sizeof(this->width));
        written_bytes += sizeof(this->width);
        out.write((char *) &this->height, sizeof(this->height));
        written_bytes += sizeof(this->height);
        out.write((char *) &this->is_extended, sizeof(this->is_extended));
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
     * @brief function to load the matching statistics supported RLPBWT object
     * @param in std::istream object from which load the matching statistics
     * supported RLPBWT structure object
     */
    void load(std::istream &in) {
        in.read((char *) &this->width, sizeof(this->width));
        in.read((char *) &this->height, sizeof(this->height));
        in.read((char *) &this->is_extended, sizeof(this->is_extended));
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

#endif  // RLPBWT_RLPBWT_NAIVE_MS_H
