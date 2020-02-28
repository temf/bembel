// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_UTIL_GENERICMATRIX_H_
#define BEMBEL_UTIL_GENERICMATRIX_H_

#include <vector>

namespace Bembel {
template <typename T>
class GenericMatrix {
 public:
  typedef typename std::vector<std::vector<T> >::size_type colIndex;
  typedef typename std::vector<T>::size_type rowIndex;
  //////////////////////////////////////////////////////////////////////////////
  //  constructors
  //////////////////////////////////////////////////////////////////////////////
  GenericMatrix() : rows_(0), cols_(0){};

  GenericMatrix(rowIndex rows, colIndex cols) { resize(rows, cols); }

  GenericMatrix(const GenericMatrix& other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    m_data_ = other.m_data_;
  }

  GenericMatrix(GenericMatrix&& other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    m_data_ = std::move(other.m_data_);
  }
  //////////////////////////////////////////////////////////////////////////////
  //  methods
  //////////////////////////////////////////////////////////////////////////////
  void resize(rowIndex rows, colIndex cols) {
    m_data_.resize(cols);
    for (auto it = m_data_.begin(); it != m_data_.end(); ++it) it->resize(rows);
    rows_ = rows;
    cols_ = cols;
    return;
  }

  colIndex cols() const { return cols_; }

  rowIndex rows() const { return rows_; }
  //////////////////////////////////////////////////////////////////////////////
  //  operators
  //////////////////////////////////////////////////////////////////////////////
  const T& operator()(rowIndex row, colIndex col) const {
    return m_data_[col][row];
  }

  T& operator()(rowIndex row, colIndex col) { return m_data_[col][row]; }

  GenericMatrix& operator=(GenericMatrix other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    std::swap(m_data_, other.m_data_);
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //  private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::vector<std::vector<T> > m_data_;
  colIndex cols_;
  rowIndex rows_;
};
}  // namespace Bembel
#endif
