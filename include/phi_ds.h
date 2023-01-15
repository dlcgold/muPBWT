//
// Created by dlcgold on 17/02/22.
//

#ifndef RLPBWT_PHI_SUPPORT_H
#define RLPBWT_PHI_SUPPORT_H


#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <optional>

#include "rl_column.h"


/**
 * @brief class to rapresent the additional data structure for phi and phi_inv
 * support
 */
class phi_ds {
private:
    /**
     * @brief default value (used to return null from the two functions)
     */
    unsigned int def{};
    /**
     * @brief default value (used to return null from the two functions)
     */
    unsigned int w{};
public:
    /**
     * @brief panel of sparse bitvectors for phi function
     */
    std::vector <sdsl::sd_vector<>> phi_vec;

    /**
     * @brief panel of sparse bitvectors for phi_inv function
     */
    std::vector <sdsl::sd_vector<>> phi_inv_vec;

    /**
     * @brief panel of rank support for phi function
     */
    std::vector <sdsl::sd_vector<>::rank_1_type> phi_rank;

    /**
     * @brief panel of rank support for phi_inv function
     */
    std::vector <sdsl::sd_vector<>::rank_1_type> phi_inv_rank;

    /**
     * @brief panel of select support for phi function
     */
    std::vector <sdsl::sd_vector<>::select_1_type> phi_select;

    /**
     * @brief panel of select support for phi_inv function
     */
    std::vector <sdsl::sd_vector<>::select_1_type> phi_inv_select;

    /**
     * @brief compressed int vector for prefix samples used by phi function
     */
    std::vector <sdsl::int_vector<>> phi_supp;

    /**
     * @brief compressed int vector for prefix samples used by phi_inv function
     */
    std::vector <sdsl::int_vector<>> phi_inv_supp;

    /**
    * @brief compressed int vector for prefix samples used by phi function
    */
    std::vector <sdsl::int_vector<>> phi_supp_l;


    /**
     * @brief default constructor
     */
    phi_ds() = default;

    /**
     * @brief default destructor
     */
    virtual ~phi_ds() = default;

    /**
     * @brief constructor of the phi/phi_inv support data structure
     * @param cols vector of the columns of the RLPBWT
     * @param panelbv random access data structure for the panel of RLPBWT
     * @param last_pref last prefix array of the PBWT
     * @param verbose bool for extra prints
     */
    explicit phi_ds(std::vector <rl_column> &cols, unsigned int h,
                    unsigned int w,
                    sdsl::int_vector<> &last_pref,
                    sdsl::int_vector<> &last_div,
                    std::vector <intv> &supp_b,
                    std::vector <intv> &supp_e,
                    bool verbose = false) {
        // default value is the panel height
        this->def = h;
        this->w = w;
        // initialize temporary panel of no-sparse bitvectors
//        auto phi_tmp = std::vector<sdsl::bit_vector>(h,
//                                                     sdsl::bit_vector(
//                                                             w + 1,
//                                                             0));
//        auto phi_inv_tmp = std::vector<sdsl::bit_vector>(h,
//                                                         sdsl::bit_vector(
//                                                                 w + 1,
//                                                                 0));
        //auto phi_tmp = std::vector<sdsl::sd_vector_builder>(h);
        //auto phi_inv_tmp = std::vector<sdsl::sd_vector_builder>(h);
        // initialize panels and vectors of the data structure
        this->phi_vec = std::vector<sdsl::sd_vector<>>(h);
        this->phi_inv_vec = std::vector<sdsl::sd_vector<>>(h);
        this->phi_rank = std::vector<sdsl::sd_vector<>::rank_1_type>(
                h);
        this->phi_inv_rank = std::vector<sdsl::sd_vector<>::rank_1_type>(
                h);
        this->phi_select = std::vector<sdsl::sd_vector<>::select_1_type>(
                h);
        this->phi_inv_select = std::vector<sdsl::sd_vector<>::select_1_type>(
                h);
        this->phi_supp = std::vector<sdsl::int_vector<>>(h);
        this->phi_inv_supp = std::vector<sdsl::int_vector<>>(h);
        this->phi_supp_l = std::vector<sdsl::int_vector<>>(h);

        // temporary vector for supports
        std::vector <std::vector<unsigned int>> phi_supp_tmp(h);
        std::vector <std::vector<unsigned int>> phi_inv_supp_tmp(h);
        std::vector <std::vector<unsigned int>> phi_supp_tmp_l(h);

        // build sparse bitvector
        for (unsigned int j = 0; j < def; j++) {
            sdsl::sd_vector_builder tmp_b = sdsl::sd_vector_builder(w + 1,
                                                                    supp_b[j].v.size());
            sdsl::sd_vector_builder tmp_e = sdsl::sd_vector_builder(w + 1,
                                                                    supp_e[j].v.size());
            for (unsigned int k = 0; k < supp_b[j].v.size(); k++) {
                tmp_b.set(supp_b[j].v[k]);
            }
            for (unsigned int k = 0; k < supp_e[j].v.size(); k++) {
                tmp_e.set(supp_e[j].v[k]);
            }
            phi_vec[j] = sdsl::sd_vector<>(tmp_b);
            phi_inv_vec[j] = sdsl::sd_vector<>(tmp_e);
        }
        for (unsigned int i = 0; i < cols.size(); i++) {
            for (unsigned int j = 0; j < cols[i].sample_beg.size(); j++) {
                // use sample beg to compute phi panel
                // use sample_end to compute
                // support phi panel (if we are in the first run we use default
                // value)
                if (j == 0) {
                    phi_supp_tmp[cols[i].sample_beg[j]].push_back(h);
                    phi_supp_tmp_l[cols[i].sample_beg[j]].push_back(0);
                } else {
                    phi_supp_tmp[cols[i].sample_beg[j]].push_back(
                            cols[i].sample_end[j - 1]);
                    phi_supp_tmp_l[cols[i].sample_beg[j]].push_back(
                            cols[i].sample_beg_lcp[j]);
                }
                // use sample end to compute phi_inv panel
                // use sample_beg to compute
                // support phi panel (if we are in the last run we use default
                // value)
                if (j == cols[i].sample_beg.size() - 1) {
                    phi_inv_supp_tmp[cols[i].sample_end[j]].push_back(
                            h);
                } else {
                    phi_inv_supp_tmp[cols[i].sample_end[j]].push_back(
                            cols[i].sample_beg[j + 1]);
                }
            }
        }
        // use the last prefix array to compute the remain values for the
        // phi support data structure (with the same "rules" of the previous
        // case)
        for (unsigned int j = 0; j < phi_supp_tmp.size(); j++) {
            if (j == 0) {
                if (phi_supp_tmp[last_pref[j]].empty() ||
                    phi_supp_tmp[last_pref[j]].back() != h) {
                    phi_supp_tmp[last_pref[j]].push_back(h);
                    phi_supp_tmp_l[last_pref[j]].push_back(0);
                }
            } else {
                if (phi_supp_tmp[last_pref[j]].empty() ||
                    phi_supp_tmp[last_pref[j]].back() != last_pref[j - 1]) {
                    phi_supp_tmp[last_pref[j]].push_back(last_pref[j - 1]);
                    phi_supp_tmp_l[last_pref[j]].push_back(last_div[j]);
                }
            }
            if (j == phi_supp_tmp.size() - 1) {
                if (phi_inv_supp_tmp[last_pref[j]].empty() ||
                    phi_inv_supp_tmp[last_pref[j]].back() != h) {
                    phi_inv_supp_tmp[last_pref[j]].push_back(h);
                }
            } else {
                if (phi_inv_supp_tmp[last_pref[j]].empty() ||
                    phi_inv_supp_tmp[last_pref[j]].back() != last_pref[j + 1]) {
                    phi_inv_supp_tmp[last_pref[j]].push_back(last_pref[j + 1]);
                }
            }
        }
        // compress and push the support sdsl int_vectors
        for (unsigned int i = 0; i < phi_supp_tmp.size(); i++) {
            sdsl::int_vector<> tmp(phi_supp_tmp[i].size());
            for (unsigned int j = 0; j < phi_supp_tmp[i].size(); j++) {
                tmp[j] = phi_supp_tmp[i][j];
            }
            sdsl::util::bit_compress(tmp);
            this->phi_supp[i] = tmp;

            sdsl::int_vector<> tmp_inv(phi_inv_supp_tmp[i].size());
            for (unsigned int j = 0; j < phi_inv_supp_tmp[i].size(); j++) {
                tmp_inv[j] = phi_inv_supp_tmp[i][j];
            }
            sdsl::util::bit_compress(tmp_inv);
            this->phi_inv_supp[i] = tmp_inv;

            sdsl::int_vector<> tmp_l(phi_supp_tmp_l[i].size());
            for (unsigned int j = 0; j < phi_supp_tmp_l[i].size(); j++) {
                tmp_l[j] = phi_supp_tmp_l[i][j];
            }
            sdsl::util::bit_compress(tmp_l);
            this->phi_supp_l[i] = tmp_l;

        }
        // create sparse bit vector and relative rank/select for phi panel
        for (unsigned int i = 0; i < def; i++) {
            this->phi_rank[i] = sdsl::sd_vector<>::rank_1_type(
                    &this->phi_vec[i]);
            this->phi_select[i] = sdsl::sd_vector<>::select_1_type(
                    &this->phi_vec[i]);
        }
        // create sparse bit vector and relative rank/select for phi_inv panel
        for (unsigned int i = 0; i < def; i++) {
            this->phi_inv_rank[i] = sdsl::sd_vector<>::rank_1_type(
                    &this->phi_inv_vec[i]);
            this->phi_inv_select[i] = sdsl::sd_vector<>::select_1_type(
                    &this->phi_inv_vec[i]);
        }


        if (verbose) {
            unsigned int index = 0;
            for (auto &col: cols) {
                std::cout << index << ":\t" << col.sample_beg
                          << "\n";
                std::cout << index << ":\t" << col.sample_end
                          << "\n";
                std::cout << "\n";
                index++;
            }
            std::cout << "\n\n";
            index = 0;
            for (auto &i: this->phi_vec) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "----------\n";
            index = 0;
            for (auto &i: this->phi_inv_vec) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "\n\n";
            index = 0;
            for (auto &i: this->phi_supp) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "----------\n";
            index = 0;
            for (auto &i: this->phi_inv_supp) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
        }
    }

    /**
     * @brief phi function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return previous prefix array value at current column (if exists)
     */
    std::optional<unsigned int> phi(unsigned int pref, unsigned int col) {
        auto res = static_cast<unsigned int>(this->phi_supp[pref][this->phi_rank[pref](
                col)]);
        if (res == this->def) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    /**
     * @brief phi_inv function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return next prefix array value at current column (if exists)
    */
    std::optional<unsigned int> phi_inv(unsigned int pref, unsigned int col) {
        auto res = static_cast<unsigned int>(this->phi_inv_supp[pref][this->phi_inv_rank[pref](
                col)]);
        if (res == this->def) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    unsigned int plcp(unsigned int pref, unsigned int col) {
        if (col == 0 || !phi(pref, col).has_value()) {
            return 0;
        }
//        unsigned int end_col = 1;
//        std::cout << this->phi_vec[pref] << std::endl;
//        std::cout << this->phi_rank[pref](col) << std::endl;
//        std::cout << this->phi_select[pref](this->phi_rank[pref](col)+1) << std::endl;
//        //std::cout << "make select\n";
//        if (col != 1) {
//            end_col = this->phi_select[pref](col - 1);
//        }

        auto end_col = this->phi_select[pref](this->phi_rank[pref](col) + 1);
        //std::cout << "end col: " << end_col << "\n";
        auto tmp = static_cast<int>(this->phi_supp_l[pref][this->phi_rank[pref](
                col)]);
        //std::cout << "value at end: " << tmp << "\n";
        if (!phi(pref, end_col + 1).has_value()) {
            //std::cout << "here\n";
            //tmp++;
        }
        if (end_col == w) {

        }
        //std::cout << "value at end: " << tmp << "\n";
        //std::cout << "diff: " << end_col - col << "\n";
        auto plcp = tmp - (end_col - col);
        return plcp;
    }

    /**
     * @brief function to obtain size in bytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false) {
        unsigned long long size = 0;
        for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
            size += sdsl::size_in_bytes(phi_vec[i]);
            size += sdsl::size_in_bytes(phi_inv_vec[i]);
            size += sdsl::size_in_bytes(phi_rank[i]);
            size += sdsl::size_in_bytes(phi_select[i]);
            size += sdsl::size_in_bytes(phi_inv_rank[i]);
            size += sdsl::size_in_bytes(phi_inv_select[i]);
            size += sdsl::size_in_bytes(phi_supp[i]);
            size += sdsl::size_in_bytes(phi_inv_supp[i]);
            size += sdsl::size_in_bytes(phi_supp_l[i]);
        }
        if (verbose) {
            std::cout << "phi support: " << size << " bytes\n";
        }
        return size;
    }

    /**
     * @brief function to obtain size in megabytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes(bool verbose = false) {
        double size_panels = 0;
        double size_supp = 0;
        for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
            size_panels += sdsl::size_in_mega_bytes(phi_vec[i]);
            size_panels += sdsl::size_in_mega_bytes(phi_inv_vec[i]);
            size_panels += sdsl::size_in_mega_bytes(phi_rank[i]);
            size_panels += sdsl::size_in_mega_bytes(phi_select[i]);
            size_panels += sdsl::size_in_mega_bytes(phi_inv_rank[i]);
            size_panels += sdsl::size_in_mega_bytes(phi_inv_select[i]);
            size_supp += sdsl::size_in_mega_bytes(phi_supp[i]);
            size_supp += sdsl::size_in_mega_bytes(phi_inv_supp[i]);
            size_supp += sdsl::size_in_mega_bytes(phi_supp_l[i]);
        }
        double size = size_panels + size_supp;
        if (verbose) {
            std::cout << "phi panels: " << size_panels << " megabytes\n";
            std::cout << "phi support: " << size_supp << " megabytes\n";
            std::cout << "phi data structure (panels + support): " << size << " megabytes\n";
        }
        return size;
    }

    /**
     * @brief function to serialize the phi/phi_inv data structure object
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "") {
        sdsl::structure_tree_node *child =
                sdsl::structure_tree::add_child(v, name,
                                                sdsl::util::class_name(
                                                        *this));
        size_t written_bytes = 0;
        out.write((char *) &this->def, sizeof(this->def));
        written_bytes += sizeof(this->def);
        out.write((char *) &this->w, sizeof(this->w));
        written_bytes += sizeof(this->w);

        for (unsigned int i = 0; i < this->phi_vec.size(); i++) {
            std::string label = "phi_vec_" + std::to_string(i);
            written_bytes += this->phi_vec[i].serialize(out, child, label);
        }

        for (unsigned int i = 0; i < this->phi_inv_vec.size(); i++) {
            std::string label = "phi_inv_vec_" + std::to_string(i);
            written_bytes += this->phi_inv_vec[i].serialize(out, child, label);
        }

        for (unsigned int i = 0; i < this->phi_supp.size(); i++) {
            std::string label = "phi_supp_" + std::to_string(i);
            written_bytes += this->phi_supp[i].serialize(out, child, label);
        }

        for (unsigned int i = 0; i < this->phi_inv_supp.size(); i++) {
            std::string label = "phi_inv_supp_" + std::to_string(i);
            written_bytes += this->phi_inv_supp[i].serialize(out, child,
                                                             label);
        }

        for (unsigned int i = 0; i < this->phi_supp_l.size(); i++) {
            std::string label = "phi_supp_l_" + std::to_string(i);
            written_bytes += this->phi_supp_l[i].serialize(out, child, label);
        }
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /**
     * @brief function to load the phi/phi_inv data structure object
     * @param in std::istream object from which load the phi/phi_inv data
     * structure object
     */
    void load(std::istream &in) {
        in.read((char *) &this->def, sizeof(this->def));
        in.read((char *) &this->w, sizeof(this->w));

        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::sd_vector<>();
            s->load(in);
            this->phi_vec.emplace_back(*s);
            delete s;
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::sd_vector<>();
            s->load(in);
            this->phi_inv_vec.emplace_back(*s);
            delete s;
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::int_vector<>();
            s->load(in);
            this->phi_supp.emplace_back(*s);
            delete s;
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::int_vector<>();
            s->load(in);
            this->phi_inv_supp.emplace_back(*s);
            delete s;
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::int_vector<>();
            s->load(in);
            this->phi_supp_l.emplace_back(*s);
            delete s;
        }
        for (auto &i: this->phi_vec) {
            this->phi_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
                    &i));
            this->phi_select.emplace_back(sdsl::sd_vector<>::select_1_type(
                    &i));
        }
        for (auto &i: this->phi_inv_vec) {
            this->phi_inv_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
                    &i));
            this->phi_inv_select.emplace_back(sdsl::sd_vector<>::select_1_type(
                    &i));
        }
    }
};


#endif //RLPBWT_PHI_SUPPORT_H
