#include <Bembel/AnsatzSpace>
#include <Bembel/H2Matrix>
#include <Bembel/IO>
#include <Bembel/Laplace>

int main() {
  Bembel::IO::Stopwatch myT;
  int refinement_level = 0;
  int polynomial_degree = 0;
  Bembel::Geometry geometry("sphere.dat");
  std::cout << "The geometry has " << geometry.get_number_of_patches()
            << " patches." << std::endl;
  myT.tic();
  Bembel::ClusterTree mesh(geometry, refinement_level);
  std::cout << "ClusterTree " << myT.toc() << std::endl;
  auto et = mesh.get_element_tree();
  Bembel::AnsatzSpace<Bembel::LaplaceSingleLayerOperator> ansatz_space(
      geometry, refinement_level, polynomial_degree);
  Bembel::LaplaceSingleLayerOperator S;
  myT.tic();
  Bembel::BlockClusterTree<double> bct0(S, ansatz_space);
  std::cout << "tree init " << myT.toc() << std::endl;

  myT.tic();
  Bembel::BlockClusterTree<double> bct =
      Bembel::BlockClusterTree<double>(S, ansatz_space);
  std::cout << "move assignment " << myT.toc() << std::endl;
  myT.tic();

  Bembel::BlockClusterTree<double> bct2 = bct0;
  std::cout << "copy assignment " << myT.toc() << std::endl;

  auto i = 0;
  std::cout << "tree copied\n";
  auto it2 = bct2.lbegin();
  for (auto it = bct.lbegin(); it != bct.lend(); ++it, ++it2) {
    assert((*it2)->get_cluster1() == (*it)->get_cluster1() &&
           (*it2)->get_cluster2() == (*it)->get_cluster2() &&
           "these should be identical!");
  }
  myT.tic();
  return 0;
}
