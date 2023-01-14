//
// Created by dlcgold on 19/02/22.
//

#include "../include/ms.h"

ms::ms() = default;
ms::~ms() = default;

ms::ms(std::vector<unsigned int> pos,
       std::vector<unsigned int> len)
        : row(std::move(pos)), len(std::move(len)) {}

ms::ms(unsigned int size) : row(size, 0), len(size, 0) {}

std::ostream &operator<<(std::ostream &os, const ms &ms) {
    os << "\nind:\t";
    for (unsigned int i = 0; i < ms.row.size(); i++) {
        os << i << "\t";
    }
    os << "\npos:\t";
    for (auto e: ms.row) {
        os << e << "\t";
    }
    os << "\nlen:\t";
    for (auto e: ms.len) {
        os << e << "\t";
    }
    return os;
}

unsigned long long ms::size_in_bytes() const {
    unsigned long long size = 0;
    size += (sizeof(unsigned int) * this->row.size() * 2);
    return size;
}

double ms::size_in_mega_bytes() const {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double)(sizeof(unsigned int) * this->row.size() * 2);
    return size * to_mega;
}
