//
// Created by dlcgold on 29/10/21.
//

#include "../include/utils.h"

char get_next_char(bool zero_first, unsigned int index_run) {
  if (zero_first) {
    if (index_run % 2 == 0) {
      return '0';
    } else {
      return '1';
    }
  } else {
    if (index_run % 2 == 0) {
      return '1';
    } else {
      return '0';
    }
  }
}

constexpr uint8_t bit_size(unsigned int value) {
  return static_cast<uint8_t>(std::floor(std::log(value)) + 1);
}

unsigned int end_interval(unsigned int start, unsigned int length) {
  return start + length - 1;
}

unsigned int u_abs(unsigned int u_a, unsigned int u_b) {
  if (u_a > u_b) {
    return u_a - u_b;
  } else {
    return u_b - u_a;
  }
}
