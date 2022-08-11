// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SRC_IO_PRINT2FILE_HPP_
#define BEMBEL_SRC_IO_PRINT2FILE_HPP_

namespace Bembel {
namespace IO {
/**
 *  \brief write Eigen::Matrix into an ascii txt file.
 **/
template <typename Derived>
int print2ascii(const std::string &fileName,
                const Eigen::MatrixBase<Derived> &var) {
  // evaluate Eigen expression into a matrix
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp =
      var;

  std::ofstream myfile;

  myfile.open(fileName);
  // write matrix to file (precision is fixed here!)
  for (auto i = 0; i < tmp.rows(); ++i) {
    for (auto j = 0; j < tmp.cols(); ++j)
      myfile << std::setprecision(10) << tmp(i, j) << " \t ";
    myfile << std::endl;
  }
  myfile.close();

  return 0;
}

/**
 *  \brief write Eigen::Matrix into a Matlab .m file.
 **/
template <typename Derived>
int print2m(const std::string &fileName, const std::string &varName,
            const Eigen::MatrixBase<Derived> &var,
            const std::string &writeMode) {
  // evaluate Eigen expression into a matrix
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp =
      var;

  std::ofstream myfile;
  // if flag is set to w, a new file is created, otherwise the new matrix
  // is just appended
  if (writeMode == "w")
    myfile.open(fileName);
  else if (writeMode == "a")
    myfile.open(fileName, std::ios_base::app);
  else
    return 1;

  myfile << varName << "=[" << std::endl;

  for (int i = 0; i < (int)tmp.rows(); ++i) {
    for (int j = 0; j < (int)tmp.cols(); ++j)
      myfile << std::setprecision(30) << tmp(i, j) << " \t ";
    myfile << std::endl;
  }
  myfile << "];" << std::endl;

  myfile.close();

  return 0;
}

template <typename Scalar>
int print2m(const std::string &fileName, const std::string &varName,
            const Eigen::SparseMatrix<Scalar> &var,
            const std::string &writeMode) {
  Eigen::VectorXd rowInd(var.nonZeros());
  Eigen::VectorXd colInd(var.nonZeros());
  Eigen::VectorXd value(var.nonZeros());
  unsigned int j = 0;
  for (auto i = 0; i < var.outerSize(); i++)
    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(var, i); it;
         ++it) {
      rowInd(j) = it.row() + 1;
      colInd(j) = it.col() + 1;
      value(j) = it.value();
      ++j;
    }
  print2m(fileName, "rows_" + varName, rowInd, writeMode);
  print2m(fileName, "cols_" + varName, colInd, "a");
  print2m(fileName, "values_" + varName, value, "a");
  std::ofstream myfile;
  // if flag is set to w, a new file is created, otherwise the new matrix
  // is just appended
  myfile.open(fileName, std::ios_base::app);
  myfile << varName << " = sparse("
         << "rows_" + varName << ","
         << "cols_" + varName << ","
         << "values_" + varName << ");\n";
  myfile.close();

  return 0;
}

template <typename Derived>
int print2bin(const std::string &fileName,
              const Eigen::MatrixBase<Derived> &var) {
  std::ofstream myfile;
  int rows = 0;
  int cols = 0;
  int IsRowMajor = 0;
  int dataSize = 0;

  IsRowMajor = var.IsRowMajor;
  dataSize = sizeof(typename Derived::Scalar);

  myfile.open(fileName, std::ios::out | std::ios::binary | std::ios::trunc);
  // write data size and row major flag
  myfile.write(reinterpret_cast<const char *>(&dataSize), sizeof(int));
  myfile.write(reinterpret_cast<const char *>(&IsRowMajor), sizeof(int));
  rows = var.rows();
  cols = var.cols();
  // write rows and cols of the matrix
  myfile.write((const char *)&(rows), sizeof(int));
  myfile.write((const char *)&(cols), sizeof(int));
  std::cout << "writing binary file" << std::endl;
  std::cout << "rows: " << rows << " cols: " << cols
            << " dataSize: " << dataSize << " IsRowMajor: " << IsRowMajor
            << std::endl;
  if (IsRowMajor) {
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic,
                  Eigen::RowMajor>
        tmp = var;

    myfile.write(reinterpret_cast<const char *>(tmp.data()),
                 rows * cols * dataSize);

  } else {
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic,
                  Eigen::ColMajor>
        tmp = var;

    myfile.write(reinterpret_cast<const char *>(tmp.data()),
                 rows * cols * dataSize);
  }

  myfile.close();

  return 0;
}

template <typename Derived>
int bin2Mat(const std::string &fileName,
            Eigen::MatrixBase<Derived> *targetMat) {
  std::ifstream myfile;
  int rows = 0;
  int cols = 0;
  int IsRowMajor = 0;
  int dataSize = 0;
  typename Derived::Scalar *data = NULL;
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      returnMat;

  myfile.open(fileName, std::ios::in | std::ios::binary);

  myfile.read(reinterpret_cast<char *>(&dataSize), sizeof(int));
  myfile.read(reinterpret_cast<char *>(&IsRowMajor), sizeof(int));
  myfile.read(reinterpret_cast<char *>(&rows), sizeof(int));
  myfile.read(reinterpret_cast<char *>(&cols), sizeof(int));
  std::cout << "reading binary file" << std::endl;
  std::cout << "rows: " << rows << " cols: " << cols
            << " dataSize: " << dataSize << " IsRowMajor: " << IsRowMajor
            << std::endl;
  if (dataSize != sizeof(typename Derived::Scalar)) {
    std::cout << "mismatch in data size of target and input file size"
              << std::endl;
    return 1;
  }

  data = new typename Derived::Scalar[rows * cols];

  myfile.read((char *)data, rows * cols * dataSize);

  myfile.close();

  if (IsRowMajor)
    returnMat =
        Eigen::Map<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                                 Eigen::Dynamic, Eigen::RowMajor> >(data, rows,
                                                                    cols);
  else
    returnMat =
        Eigen::Map<Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                                 Eigen::Dynamic, Eigen::ColMajor> >(data, rows,
                                                                    cols);

  *targetMat = returnMat;

  delete[] data;
  return 0;
}
}  // namespace IO
}  // namespace Bembel

#endif  // BEMBEL_SRC_IO_PRINT2FILE_HPP_
