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
