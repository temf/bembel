#ifndef __WRITEMESH2VTK__
#define __WRITEMESH2VTK__

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

/** \brief
 *
 */
namespace Eigen {
template <typename Derived1, typename Derived2, typename Derived3>
void writeMesh2vtk(const std::string &fileName,
                   const Eigen::MatrixBase<Derived1> &P,
                   const Eigen::MatrixBase<Derived2> &E,
                   const Eigen::MatrixBase<Derived3> &Cdata) {
  std::ofstream myfile;
  myfile.open(fileName);
  myfile << "# vtk DataFile Version 3.1\n";
  myfile << "this file hopefully represents my surface now\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";
  // print point list
  myfile << "POINTS " << P.cols() << " FLOAT\n";
  for (auto i = 0; i < P.cols(); ++i)
    myfile << float(P(0, i)) << " " << float(P(1, i)) << " " << float(P(2, i))
           << "\n";
  myfile << "\n";

  // print element list
  myfile << "CELLS " << E.cols() << " " << 5 * E.cols() << "\n";
  for (auto i = 0; i < E.cols(); ++i)
    myfile << int(4) << " " << int(E(0, i)) << " " << int(E(1, i)) << " "
           << int(E(2, i)) << " " << int(E(3, i)) << "\n";
  myfile << "\n";

  myfile << "CELL_TYPES " << E.cols() << "\n";
  for (auto i = 0; i < E.cols(); ++i) myfile << int(9) << "\n";
  myfile << "\n";
  // print cell labels
  myfile << "POINT_DATA " << P.cols() << "\n";
  myfile << "SCALARS value FLOAT\n";
  myfile << "LOOKUP_TABLE default\n";
  for (auto i = 0; i < P.cols(); ++i) myfile << float(Cdata(i)) << "\n";

  myfile.close();
#if 0
/* print z-values of the geometry and solved density for visualization        */
fprintf (f, "POINT_DATA %d\n", np);
fprintf (f, "SCALARS Solution FLOAT\n");
fprintf (f, "LOOKUP_TABLE default\n");
for (i = 0 ; i < np ; ++i)
   fprintf (f, "%20.16f\n", u[i]);
fprintf (f, "\n");
#endif
  return;
}
}  // namespace Eigen
#endif
