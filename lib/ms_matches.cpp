//
// Created by dlcgold on 20/02/22.
//

#include "../include/ms_matches.h"

ms_matches::ms_matches() = default;


std::ostream &operator<<(std::ostream &os, const ms_matches &ms) {
    os << "\nmatches (<ending>, <length>, [haplotypes]):\n";
    bool haplo = ms.haplos.empty();
    for (unsigned int i = 0; i < ms.basic_matches.size(); i++) {
        os << std::get<2>(ms.basic_matches[i]) << ", " << std::get<1>(ms.basic_matches[i]);
        if (!haplo) {
            os << ", [";
            for (auto h: ms.haplos[i]) {
                os << h << " ";
            }
            os << "]";
        }
        os << "\n";
    }
    return os;
}



unsigned long long ms_matches::size_in_bytes() const {
    unsigned long long size = 0;
    bool haplo = this->haplos.empty();
    size += (sizeof(unsigned int) * this->basic_matches.size() * 3);
    if(!haplo){
        unsigned int nhaplo = 0;
        for(const auto& h: haplos){
            nhaplo += h.size();
        }
        size += (sizeof(unsigned int) * nhaplo);
    }
    return size;
}

double ms_matches::size_in_mega_bytes() const {
    double size = 0;
    bool haplo = haplos.empty();
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double)(sizeof(unsigned int) * this->basic_matches.size() * 3);
    if(!haplo){
        unsigned int nhaplo = 0;
        for(const auto& h: haplos){
            nhaplo += h.size();
        }
        size += (double)(sizeof(unsigned int) * nhaplo);
    }
    return size * to_mega;
}