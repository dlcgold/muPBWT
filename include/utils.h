//
// Created by dlcgold on 29/10/21.
//

#ifndef RLPBWT_UTILS_H
#define RLPBWT_UTILS_H

#include <climits>
#include <sdsl/config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>
#include <string>
#include <vector>

class intv {
public:
  int size;
  int capacity;
  // sdsl::int_vector<t_width> v;
  sdsl::int_vector<> v;
  intv() {
    size = 0;
    capacity = 0;
    //        v.resize(1);
    //        v[0]=max_value;
    //        sdsl::util::bit_compress(v);
    //        v.resize(0);
  }

  void push_back(const unsigned long int &x) {
    if (size == capacity) {
      capacity = size * 2 + 1;
      v.resize(capacity);
    }
    v[size++] = x;
  }

  void compress() {
    v.resize(size);
    sdsl::util::bit_compress(v);
  }
};

/**
 * @brief function to get the char (in bialleic case with 0 and 1) at certain
 * run index
 * @param zero_first bool to check first value of the column
 * @param index_run run index in the run-length encoded PBWT column
 * @return the char at the queried run
 */
char get_next_char(bool zero_first, unsigned int index_run);
int get_intervall(unsigned int start, unsigned int size);
constexpr uint8_t bit_size(unsigned int value);

template <typename T> double vectorsizeof(const typename std::vector<T> &vec) {
  return (sizeof(T) * vec.size()) * 0.000001;
}

inline unsigned int dense_size_byte(unsigned int h, unsigned int w) {
  sdsl::int_vector<> ext_size(1, h);
  unsigned int size = 0;
  size += (sdsl::size_in_bytes(ext_size) * h * w);
  size = (size * 5) + (sdsl::size_in_bytes(ext_size) * w);
  return size;
}

inline double dense_size_megabyte(unsigned int h, unsigned int w) {
  sdsl::int_vector<> ext_size(1, h);
  double size = 0;
  size += (sdsl::size_in_mega_bytes(ext_size) * h * w);
  size = (size * 5) + (sdsl::size_in_mega_bytes(ext_size) * w);
  return size;
}

// THANKS to Massimiliano Rossi:
// https://github.com/maxrossi91/moni/blob/595da8cb01376074ba74e13273fc9072f5af410f/include/common/common.hpp
template <class T, typename size_type>
uint64_t my_serialize_array(
    const T *p, const size_type size, std::ostream &out,
    typename std::enable_if<std::is_fundamental<T>::value>::type * = 0) {
  size_t written_bytes = 0;
  if (size > 0) {
    size_type idx = 0;
    while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size)) {
      out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
      written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
      p += sdsl::conf::SDSL_BLOCK_SIZE;
      idx += sdsl::conf::SDSL_BLOCK_SIZE;
    }
    out.write((char *)p, ((size)-idx) * sizeof(T));
    written_bytes += ((size)-idx) * sizeof(T);
  }
  return written_bytes;
}

//! Serialize each element of an std::vector
/*!
 * \param vec The vector which should be serialized.
 * \param out Output stream to which should be written.
 * \param v   Structure tree node. Note: If all elements have the same
 *            structure, then it is tried to combine all elements (i.e.
 *            make one node w with size set to the cumulative sum of all
 *           sizes of the children)
 */
// specialization for fundamental types
template <class T>
uint64_t my_serialize_vector(
    const std::vector<T> &vec, std::ostream &out, sdsl::structure_tree_node *v,
    const std::string &name,
    typename std::enable_if<std::is_fundamental<T>::value>::type * = 0) {
  if (!vec.empty()) {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(
        v, name, "std::vector<" + sdsl::util::class_name(vec[0]) + ">");
    size_t written_bytes = 0;

    const T *p = &vec[0];
    typename std::vector<T>::size_type idx = 0;
    while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size())) {
      out.write((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
      written_bytes += sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T);
      p += sdsl::conf::SDSL_BLOCK_SIZE;
      idx += sdsl::conf::SDSL_BLOCK_SIZE;
    }
    out.write((char *)p, ((vec.size()) - idx) * sizeof(T));
    written_bytes += ((vec.size()) - idx) * sizeof(T);

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  } else {
    return 0;
  }
}

/**
 * @brief function to serialize custom objects
 * @tparam X type of the object
 * @param x object to serialize
 * @param in std::ostream object in which serialize the data
 */
template <typename X>
uint64_t my_serialize(
    const std::vector<X> &x, std::ostream &out,
    sdsl::structure_tree_node *v = nullptr, std::string name = "",
    typename std::enable_if<std::is_fundamental<X>::value>::type * = 0) {
  return sdsl::serialize(x.size(), out, v, name) +
         my_serialize_vector(x, out, v, name);
}

/**
 * @brief Load an array of size elements into p. p should be preallocated.
 *
 * \tparam T
 * \tparam size_type
 * @param p
 * @param size
 * @param in
 */
template <class T, typename size_type>
void my_load_array(
    T *p, const size_type size, std::istream &in,
    typename std::enable_if<std::is_fundamental<T>::value>::type * = 0) {
  size_type idx = 0;
  while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (size)) {
    in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    p += sdsl::conf::SDSL_BLOCK_SIZE;
    idx += sdsl::conf::SDSL_BLOCK_SIZE;
  }
  in.read((char *)p, ((size)-idx) * sizeof(T));
}

//! Load all elements of a vector from a input stream
/*! \param vec  Vector whose elements should be loaded.
 *  \param in   Input stream.
 *  \par Note
 *   The vector has to be resized prior the loading
 *   of its elements.
 */
template <class T>
void my_load_vector(
    std::vector<T> &vec, std::istream &in,
    typename std::enable_if<std::is_fundamental<T>::value>::type * = 0) {
  T *p = &vec[0];
  typename std::vector<T>::size_type idx = 0;
  while (idx + sdsl::conf::SDSL_BLOCK_SIZE < (vec.size())) {
    in.read((char *)p, sdsl::conf::SDSL_BLOCK_SIZE * sizeof(T));
    p += sdsl::conf::SDSL_BLOCK_SIZE;
    idx += sdsl::conf::SDSL_BLOCK_SIZE;
  }
  in.read((char *)p, ((vec.size()) - idx) * sizeof(T));
}

/**
 * @brief function to load custom objects
 * @tparam X type of the object
 * @param x object in which load the data
 * @param in std::istream object from which load the data
 */
template <typename X>
void my_load(
    std::vector<X> &x, std::istream &in,
    typename std::enable_if<std::is_fundamental<X>::value>::type * = 0) {
  typename std::vector<X>::size_type size;
  sdsl::load(size, in);
  x.resize(size);
  my_load_vector(x, in);
}

#endif // RLPBWT_UTILS_H
