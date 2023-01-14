//
// Created by dlcgold on 30/10/21.
//

#ifndef RLPBWT_EXCEPTIONS_H
#define RLPBWT_EXCEPTIONS_H

#include <exception>

/**
 * @brief exception file not found
 */
class FileNotFoundException: public std::exception {
  const char* what() const noexcept override {
    return "file not found";
  }
};

/**
 * @brief exception file not well formed
 */
class FileNotGoodException: public std::exception {
    const char* what() const noexcept override {
        return "file not good for queries";
    }
};

/**
 * @brief exception query of different length
 */
class NotEqualLengthException: public std::exception {
    const char* what() const noexcept override {
        return "Query has not the same length of panel";
    }
};

/**
 * @brief exception file slp not found
 */
class SlpNotFoundException: public std::exception {
    const char* what() const noexcept override {
        return "SLP file is required";
    }
};

/**
 * @brief exception wrong RA type
 */
class WrongRaTypeException: public std::exception {
    const char* what() const noexcept override {
        return "rlpbwt type could be panel_ra or slp_panel_ra";
    }
};

/**
 * @brief exception RLPBWT not build with thresholds support
 */
class NoThrException: public std::exception {
    const char* what() const noexcept override {
        return "rlpbwt is not build with thresholds";
    }
};

#endif //RLPBWT_EXCEPTIONS_H
