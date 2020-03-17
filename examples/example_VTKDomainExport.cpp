#include <Bembel/IO>

int main() {
  using namespace Bembel;
  // The VTKwriter sets up initial geomety information.
  Eigen::VectorXd linspace = Eigen::VectorXd::LinSpaced(50, 0, 3);
  VTKDomainExport writer(linspace, linspace, linspace);

  Geometry geo("sphere.dat");

  // Now we can add user defined data. There are different options.
  // Exceptionally handy is a dataset that describes the distance to the
  // geometry. This is useful for paraview visualizations, s.t. one can use the
  // threashold feature to clip the domain. For this, we can just iterate over
  // the elements of a ClusterTree. Usually the ClusterTree is hidden in an ansatz space and
  // can be accest via the get_mesh() method. Try running the
  // example_VTKSurfaceExport and visualise both files on top of each other.
  ClusterTree mesh(geo, 5);
  std::function<double(const Eigen::Vector3d&)> fun1 =
      [&](const Eigen::Vector3d& point_in_space) {
        double dist = 1000;
        for (auto e = mesh.get_element_tree().cpbegin();
             e != mesh.get_element_tree().cpend(); ++e) {
          double new_dist = ((*e).midpoint_ - point_in_space).norm();
          dist = dist > new_dist ? new_dist : dist;
        }
        return dist;
      };

  // One can also visualize vector fields.

  std::function<Eigen::Vector3d(const Eigen::Vector3d&)> fun2 =
      [](const Eigen::Vector3d& point_in_space) {
        return Eigen::Vector3d(point_in_space(0),
                               .5 * point_in_space(2) * point_in_space(1),
                               point_in_space(0));
      };

  writer.addDataSet("Distance_to_sphere", fun1);
  writer.addDataSet("Some_vector_field", fun2);

  // Finally, we print to file.
  writer.writeToFile("example.vts");
  return 0;
}
