//
// Created by dlcgold on 16/08/22.
//

#include "../include/rl_column.h"

rl_column::rl_column(bool zeroFirst, unsigned int count0,
                     sdsl::int_vector<> p,
                     sdsl::int_vector<> uv,
                     sdsl::int_vector<> t,
                     sdsl::int_vector<> sample_beg,
                     sdsl::int_vector<> sample_end,
                     sdsl::int_vector<> sample_beg_lcp) :
        zero_first(zeroFirst),
        count_0(count0),
        p(std::move(p)),
        uv(std::move(uv)),
        t(std::move(t)),
        sample_beg(std::move(sample_beg)),
        sample_end(std::move(sample_end)),
        sample_beg_lcp(std::move(sample_beg_lcp)) {}


rl_column::rl_column() = default;

rl_column::~rl_column() = default;

std::ostream &operator<<(std::ostream &os, const rl_column &naive) {
    std::string value;
    if (naive.zero_first) {
        value = "0";
    } else {
        value = "1";
    }
    os << "first value: " << value << "\n";
    os << "c value: " << naive.count_0 << "\n";
    if (naive.t.empty()) {
        for (unsigned int i = 0; i < naive.p.size(); i++) {
            os << naive.p[i] << "\t" << naive.uv[i] << "\t("
               << naive.sample_beg[i] << ", " << naive.sample_end[i]
               << ")\n";
        }
    } else {
        for (unsigned int i = 0; i < naive.p.size(); i++) {
            os << naive.p[i] << "\t" << naive.uv[i] << "\t" << naive.t[i]
               << "\t(" << naive.sample_beg[i] << ", " << naive.sample_end[i]
               << ")\t" << value << "\n";
            if (value[0] == '0') {
                value[0] = '1';
            } else {
                value[0] = '0';
            }
        }

    }
    return os;
}

unsigned long long rl_column::size_in_bytes(bool verbose) const {
    unsigned long long size = 0;
    size += sizeof(bool);
    size += sizeof(unsigned int);
    auto p_size = sdsl::size_in_bytes(this->p);
    auto uv_size = sdsl::size_in_bytes(this->uv);
    auto t_size = sdsl::size_in_bytes(this->t);
    auto sample_size = sdsl::size_in_bytes(this->sample_beg) +
                       sdsl::size_in_bytes(this->sample_end) +
                       sdsl::size_in_bytes(this->sample_beg_lcp);
    size += p_size;
    size += uv_size;
    size += t_size;
    size += sample_size;
    if (verbose) {
        std::cout << "p: " << p_size << " bytes\n";
        std::cout << "uv: " << uv_size << " bytes\n";
        std::cout << "thr: " << uv_size << " bytes\n";
        std::cout << "sample: " << sample_size << " bytes\n";
    }
    return size;
}

double rl_column::size_in_mega_bytes(bool verbose) const {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double) (sizeof(bool) * to_mega);
    size += (double) (sizeof(unsigned int) * to_mega);
    auto p_size = sdsl::size_in_mega_bytes(this->p);
    auto uv_size = sdsl::size_in_mega_bytes(this->uv);
    auto t_size = sdsl::size_in_mega_bytes(this->t);
    auto sample_size = sdsl::size_in_mega_bytes(this->sample_beg) +
                       sdsl::size_in_mega_bytes(this->sample_end) +
                       sdsl::size_in_mega_bytes(this->sample_beg_lcp);
    size += p_size;
    size += uv_size;
    size += t_size;
    size += sample_size;

    if (verbose) {
        std::cout << "p: " << p_size << " megabytes\n";
        std::cout << "uv: " << uv_size << " megabytes\n";
        std::cout << "thr: " << uv_size << " bytes\n";
        std::cout << "sample: " << sample_size << " bytes\n";
    }
    return size;
}

size_t
rl_column::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                     const std::string &name) {
    sdsl::structure_tree_node *child =
            sdsl::structure_tree::add_child(v, name,
                                            sdsl::util::class_name(
                                                    *this));
    size_t written_bytes = 0;

    out.write((char *) &this->zero_first, sizeof(this->zero_first));
    written_bytes += sizeof(this->zero_first);

    out.write((char *) &this->count_0, sizeof(this->count_0));
    written_bytes += sizeof(this->count_0);

    this->p.serialize(out, child, "p");
    this->uv.serialize(out, child, "uv");
    this->t.serialize(out, child, "t");
    this->sample_beg.serialize(out, child, "s_beg");
    this->sample_end.serialize(out, child, "s_end");
    //this->sample_beg_lcp.serialize(out, child, "s_beg_l");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void rl_column::load(std::istream &in) {
    in.read((char *) &this->zero_first, sizeof(this->zero_first));
    in.read((char *) &this->count_0, sizeof(this->count_0));
    this->p.load(in);
    this->uv.load(in);
    this->t.load(in);
    this->sample_beg.load(in);
    this->sample_end.load(in);
    //this->sample_beg_lcp.load(in);
}




