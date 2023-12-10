//
// Created by dlcgold on 19/02/22.
//

#include "../include/ms_k.h"

ms_k::ms_k() = default;
ms_k::~ms_k() = default;

ms_k::ms_k(std::vector<unsigned int> row, std::vector<unsigned int> len,
           std::vector<unsigned int> row_supp,
           std::vector<unsigned int> len_supp)
    : row(std::move(row)), len(std::move(len)), row_supp(std::move(row_supp)),
      len_supp(std::move(len_supp)) {}

ms_k::ms_k(unsigned int size)
    : row(size, 0), len(size, 0), row_supp(size, 0), len_supp(size, 0) {}

std::ostream &operator<<(std::ostream &os, const ms_k &ms_k) {
  os << "\nind:\t\t";
  for (unsigned int i = 0; i < ms_k.row.size(); i++) {
    os << i << "\t";
  }
  os << "\npos:\t\t";
  for (auto e : ms_k.row) {
    os << e << "\t";
  }
  os << "\nlen:\t\t";
  for (auto e : ms_k.len) {
    os << e << "\t";
  }
  os << "\npos_supp:\t";
  for (auto e : ms_k.row_supp) {
    os << e << "\t";
  }
  os << "\nlen_supp:\t";
  for (auto e : ms_k.len_supp) {
    os << e << "\t";
  }
  return os;
}

unsigned long long ms_k::size_in_bytes() const {
  unsigned long long size = 0;
  size += (sizeof(unsigned int) * this->row.size() * 4);
  return size;
}

double ms_k::size_in_mega_bytes() const {
  double size = 0;
  double to_mega = ((double)1 / (double)1024) / (double)1024;
  size += (double)(sizeof(unsigned int) * this->row.size() * 4);
  return size * to_mega;
}
