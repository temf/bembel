// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_IO_VTKPOINTEXPORT_H_
#define BEMBEL_IO_VTKPOINTEXPORT_H_

namespace Bembel {

// This class provides the possibilty to generate a vtk-visualization.
class VTKPointExport {
 public:
  /**
   * \ingroup IO
   * \brief Provides export routines to the VTK file format.
   **/

  VTKPointExport(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    init_VTKPointExport(points);
  }

  inline void init_VTKPointExport(
      const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    points_ = points;
    point_number_ = Eigen::VectorXi(points_.rows());
    for (int i = 0; i < points_.rows(); ++i) point_number_(i) = i;
    return;
  }

  // This routine turns the data of the above DataSet-routines into a string and
  // stores it.
  inline void addDataSet(const std::string &name, const Eigen::MatrixXd &mat) {
    assert(mat.cols() == 1 || mat.cols() == 3);

    const int cols = mat.cols();
    std::string data_ascii = "<DataArray type=\"Float32\" Name=\"" + name +
                             "\" NumberOfComponents=\"" +
                             std::to_string(mat.cols()) +
                             "\" format=\"ascii\">\n";
    for (int i = 0; i < mat.rows(); ++i) {
      for (int j = 0; j < cols; ++j) {
        data_ascii.append(std::to_string(mat(i, j)) + " ");
      }
      data_ascii.append("\n");
    }
    data_ascii.append("</DataArray>\n");
    additionalData_.push_back(data_ascii);
  }

  inline void writeToFile(const std::string &filename) {
    std::ofstream output;
    output.open(filename);
    output << "<VTKFile type=\"PolyData\" version=\"0.1\" "
              "byte_order=\"LittleEndian\">\n"
              "<PolyData>\n"
              "<Piece NumberOfPoints=\" "
           << points_.rows()
           << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" "
              "NumberOfPolys=\"0\">\n"
              "<Points>\n"
              "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
              "format=\"ascii\">\n";
    output << points_;
    output << "\n</DataArray>\n</Points>\n";
    output << "<PointData Scalars=\"number\">\n"
              "<DataArray type=\"Int32\" Name=\"number\" "
              "format=\"ascii\">\n";
    output << point_number_;
    output << "</DataArray>\n";
    for (auto data : additionalData_) {
      output << data;
    }
    output
        << "</PointData>\n"
           "</Piece>\n"
           "</PolyData>\n"
           "</VTKFile>";
    output.close();
    return;
  }

  // can be used to clear the additional Data.
  inline void clearData() {
    additionalData_ = {};
    additionalData_.shrink_to_fit();
  }

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 3> points_;
  Eigen::VectorXi point_number_;
  std::vector<std::string> additionalData_;
};  // namespace Bembel

}  // namespace Bembel

#endif
