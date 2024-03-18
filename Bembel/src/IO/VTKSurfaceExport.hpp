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

#ifndef BEMBEL_SRC_IO_VTKSURFACEEXPORT_HPP_
#define BEMBEL_SRC_IO_VTKSURFACEEXPORT_HPP_

namespace Bembel {

/**
 * \ingroup IO
 * \brief Provides export routines from functions on geometries to the VTK file
 * format.
 *
 * One can add data to visualize via the addDataSet methods.
 */
class VTKSurfaceExport {
 public:
  /**
   * The constructor wants a geometetry and a refinement level. This choice is
   * deliberately not a mesh, since the visualization will often be on a finer
   * mesh then that of a computation.
   */
  VTKSurfaceExport(const Geometry& geo, int M) {
    init_VTKSurfaceExport(geo, M);
  }
  inline void init_VTKSurfaceExport(const Geometry& geo, int M) {
    msh.init_ClusterTree(geo, M);
    points = msh.get_element_tree().generatePointList().transpose();
    cells = msh.get_element_tree().generateElementList().transpose();
    normals = Eigen::MatrixXd(cells.rows(), 3);
    patch_number = Eigen::VectorXi(cells.rows());

    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      normals.row(e->id_) = msh.get_geometry()[e->patch_]
                                .evalNormal(e->referenceMidpoint())
                                .normalized()
                                .transpose();
    }
    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      patch_number(e->id_) = e->patch_;
    }
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<double(int, const Eigen::Vector2d&)> fun) {
    addDataSet_(name, getData_<double>(fun));
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<std::complex<double>(int, const Eigen::Vector2d&)> fun) {
    Eigen::VectorXcd data = getData_<std::complex<double>>(fun);
    addDataSet_(name + std::string("_real"), data.real());
    addDataSet_(name + std::string("_imag"), data.imag());
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<Eigen::Vector3d(int, const Eigen::Vector2d&)> fun) {
    addDataSet_(name, getData_<double, 3>(fun));
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<Eigen::Vector3cd(int, const Eigen::Vector2d&)> fun) {
    Eigen::MatrixXcd data = getData_<std::complex<double>, 3>(fun);
    addDataSet_(name + std::string("_real"), data.real());
    addDataSet_(name + std::string("_imag"), data.imag());
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<Eigen::Vector3d(const Eigen::Vector3d&)> fun) {
    addDataSet_(name, getData_<double, 3>(fun));
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<Eigen::Vector3cd(const Eigen::Vector3d&)> fun) {
    Eigen::MatrixXcd data = getData_<std::complex<double>, 3>(fun);
    addDataSet_(name + std::string("_real"), data.real());
    addDataSet_(name + std::string("_imag"), data.imag());
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(const std::string& name,
                         std::function<double(const Eigen::Vector3d&)> fun) {
    addDataSet_(name, getData_<double>(fun));
    return;
  }
  /**
   * \brief Add a function of the given type to the visualization.
   */
  inline void addDataSet(
      const std::string& name,
      std::function<std::complex<double>(const Eigen::Vector3d&)> fun) {
    Eigen::VectorXcd data = getData_<std::complex<double>>(fun);
    addDataSet_(name + std::string("_real"), data.real());
    addDataSet_(name + std::string("_imag"), data.imag());
    return;
  }
  /**
   * \brief Add a function given through coefficients of an AnsatzSpace to the
   * visualization.
   */
  template <typename LinOp>
  inline void addDataSet(
      const std::string& name, const AnsatzSpace<LinOp>& ansatz_space,
      const Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                          Eigen::Dynamic, Eigen::Dynamic>& coefficients) {
    // Initialize FunctionEvaluator
    FunctionEvaluator<LinOp> evaluator(ansatz_space);
    evaluator.set_function(coefficients);

    // Define function for FunctionEvaluator
    typedef typename DifferentialFormTraits<
        LinearOperatorTraits<LinOp>::Form,
        typename LinearOperatorTraits<LinOp>::Scalar>::FunctionSpaceValue
        ReturnType;
    constexpr int dimension =
        DifferentialFormTraits<LinearOperatorTraits<LinOp>::Form,
                               typename LinearOperatorTraits<LinOp>::Scalar>::
            FunctionSpaceOutputDimension;
    std::function<ReturnType(int, const Eigen::Vector2d&)> fun =
        getEvaluatorFunction_<LinOp, dimension>::get(evaluator);

    // Define function for FunctionEvaluator
    addDataSet(name, fun);

    return;
  }

  /**
   * \brief Write geometry with passed visualization data to file.
   */
  inline void writeToFile(const std::string& filename) {
    std::ofstream output;
    output.open(filename);
    output << "<VTKFile type=\"PolyData\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n"
              "<PolyData>\n"
              "<Piece NumberOfPoints=\" "
           << points.rows()
           << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" "
              "NumberOfPolys=\""
           << cells.rows()
           << "\">\n"
              "<Points>\n"
              "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    output << points;
    output << "</DataArray>\n</Points>";
    output << "<CellData Scalars=\"patch_number\" Normals=\"cell_normals\">\n"
              "<DataArray type=\"Int32\" Name=\"patch_number\" "
              "format=\"ascii\">\n";
    output << patch_number;
    output << "</DataArray>\n<DataArray type=\"Float32\" Name=\"cell_normals\" "
              "NumberOfComponents=\"3\" format=\"ascii\">";
    output << normals;
    output << "</DataArray>\n";
    for (auto data : additionalData) {
      output << data;
    }
    output << "</CellData>\n"
              "<Polys>\n"
              "<DataArray type=\"Int32\" Name=\"connectivity\" "
              "format=\"ascii\">\n";
    output << cells;
    output << "</DataArray>\n"
              "<DataArray type=\"Int32\" Name=\"offsets\" "
              "format=\"ascii\">\n";
    for (int i = 0; i < cells.rows(); ++i) {
      output << (i + 1) * 4;
      output << "\n";
    }
    output << "</DataArray>\n"
              "</Polys>\n"
              "</Piece>\n"
              "</PolyData>\n"
              "</VTKFile>";
    output.close();
    return;
  }
  // can be used to clear the additional Data.
  inline void clearData() {
    additionalData = {};
    additionalData.shrink_to_fit();
  }

 private:
  // This routine turns the data of the above DataSet-routines into a string and
  // stores it.
  inline void addDataSet_(const std::string& name, const Eigen::MatrixXd& mat) {
    assert(mat.cols() == 1 || mat.cols() == 3);
    assert(mat.rows() == cells.rows());

    const int cols = mat.cols();
    std::string data_ascii = "<DataArray type=\"Float32\" Name=\"" + name +
                             "\" NumberOfComponents=\"" +
                             std::to_string(mat.cols()) +
                             "\" format=\"ascii\">\n";
    std::ostringstream out;
    out.precision(6);
    for (int i = 0; i < mat.rows(); ++i) {
      for (int j = 0; j < cols; ++j) {
        out << std::scientific << mat(i, j);
        data_ascii.append(std::move(out).str() + " ");
        out.str("");
        out.clear();
      }
      data_ascii.append("\n");
    }
    data_ascii.append("</DataArray>\n");
    additionalData.push_back(data_ascii);
  }

  // generate data for functions living in spatial domains
  template <typename Scalar, int dim>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, dim> getData_(
      const std::function<
          Eigen::Matrix<Scalar, dim, 1>(const Eigen::Vector3d&)>& fun) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> data(cells.rows(),
                                                               dim);
    int k = 0;
    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      data.row(e->id_) =
          fun(msh.get_geometry()[e->patch_].eval(e->referenceMidpoint()))
              .transpose();
    }
    return data;
  }
  template <typename Scalar>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, 1> getData_(
      std::function<Scalar(const Eigen::Vector3d&)>& fun) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> data(cells.rows());
    int k = 0;
    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      data(e->id_) =
          fun(msh.get_geometry()[e->patch_].eval(e->referenceMidpoint()));
    }
    return data;
  }

  // generate data for functions living on patches
  template <typename Scalar, int dim>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, dim> getData_(
      const std::function<
          Eigen::Matrix<Scalar, dim, 1>(int, const Eigen::Vector2d&)>& fun) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> data(cells.rows(),
                                                               dim);
    int k = 0;
    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      data.row(e->id_) = fun(e->patch_, e->referenceMidpoint()).transpose();
    }
    return data;
  }
  template <typename Scalar>
  inline Eigen::Matrix<Scalar, Eigen::Dynamic, 1> getData_(
      std::function<Scalar(int, const Eigen::Vector2d&)>& fun) {
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> data(cells.rows());
    int k = 0;
    for (auto e = msh.get_element_tree().cpbegin();
         e != msh.get_element_tree().cpend(); ++e) {
      data(e->id_) = fun(e->patch_, e->referenceMidpoint());
    }
    return data;
  }

  // helper function for FunctionEvaluator based output
  // template <typename LinOp, typename ReturnType>
  template <typename LinOp, int size>
  struct getEvaluatorFunction_ {
    static std::function<
        Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, size, 1>(
            int, const Eigen::Vector2d&)>
    get(const FunctionEvaluator<LinOp>& evaluator) {
      typedef typename LinearOperatorTraits<LinOp>::Scalar Scalar;
      typedef
          typename Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                                 size, 1>
              ReturnType;
      std::function<ReturnType(int, const Eigen::Vector2d&)> fun =
          [&](int patch_number, const Eigen::Vector2d& reference_domain_point) {
            return evaluator.evaluateOnPatch(patch_number,
                                             reference_domain_point);
          };
      return fun;
    }
  };
  template <typename LinOp>
  struct getEvaluatorFunction_<LinOp, 1> {
    static std::function<typename LinearOperatorTraits<LinOp>::Scalar(
        int, const Eigen::Vector2d&)>
    get(const FunctionEvaluator<LinOp>& evaluator) {
      typedef typename LinearOperatorTraits<LinOp>::Scalar Scalar;
      std::function<typename DifferentialFormTraits<
          LinearOperatorTraits<LinOp>::Form,
          typename LinearOperatorTraits<LinOp>::Scalar>::
                        FunctionSpaceValue(int, const Eigen::Vector2d&)>
          fun = [&](int patch_number,
                    const Eigen::Vector2d& reference_domain_point) {
            return evaluator.evaluateOnPatch(patch_number,
                                             reference_domain_point)(0);
          };
      return fun;
    }
  };

  ClusterTree msh;
  Eigen::MatrixXd points;
  Eigen::MatrixXi cells;
  Eigen::MatrixXd normals;
  Eigen::VectorXi patch_number;
  std::vector<std::string> additionalData;
};

}  // namespace Bembel

#endif  // BEMBEL_SRC_IO_VTKSURFACEEXPORT_HPP_
