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
//
#include <Bembel/AnsatzSpace>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/Laplace>

int main() {
  Bembel::IO::Stopwatch sw;
  Bembel::Geometry geometry("sphere.dat");
  std::cout << std::string(60, '-') << std::endl;
  std::cout << "The geometry has " << geometry.get_number_of_patches()
            << " patches." << std::endl;
  std::cout << std::string(60, '-') << std::endl;
  std::cout << "    cluster tree setup\n";
  for (int refinement_level = 0; refinement_level < 11; ++refinement_level) {
    sw.tic();
    Bembel::ClusterTree mesh(geometry, refinement_level);
    std::cout << std::setw(2) << refinement_level
              << ") setup time: " << sw.toc() << std::endl;
  }
  std::cout << std::string(60, '-') << std::endl;
  std::cout << "    block cluster tree setup\n";
  for (int refinement_level = 0; refinement_level < 5; ++refinement_level) {
    sw.tic();
    const int polynomial_degree = 3;
    Bembel::AnsatzSpace<Bembel::LaplaceSingleLayerOperator> ansatz_space(
        geometry, refinement_level, polynomial_degree);
    Bembel::LaplaceSingleLayerOperator S;
    Bembel::BlockClusterTree<double> bct0(S, ansatz_space);
    std::cout << std::setw(2) << refinement_level
              << ") setup time: " << sw.toc() << std::endl;
  }
  std::cout << std::string(60, '-') << std::endl;
  {
    const int polynomial_degree = 3;
    const int refinement_level = 5;
    Bembel::AnsatzSpace<Bembel::LaplaceSingleLayerOperator> ansatz_space(
        geometry, refinement_level, polynomial_degree);
    Bembel::LaplaceSingleLayerOperator S;
    Bembel::BlockClusterTree<double> bct0(S, ansatz_space);
    sw.tic();
    Bembel::BlockClusterTree<double> bct =
        std::move(Bembel::BlockClusterTree<double>(S, ansatz_space));
    std::cout << "move assignment " << sw.toc() << std::endl;
    sw.tic();

    Bembel::BlockClusterTree<double> bct2 = bct0;
    std::cout << "copy assignment " << sw.toc() << std::endl;

    auto i = 0;
    std::cout << "tree copied\n";
    auto it2 = bct2.lbegin();
    for (auto it = bct.lbegin(); it != bct.lend(); ++it, ++it2) {
      assert((*it2)->get_cluster1() == (*it)->get_cluster1() &&
             (*it2)->get_cluster2() == (*it)->get_cluster2() &&
             "these should be identical!");
    }
  }
  return 0;
}
