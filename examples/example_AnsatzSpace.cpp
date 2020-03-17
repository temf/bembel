#include <Bembel/AnsatzSpace>
#include <Bembel/DummyOperator>
#include "Bembel/src/AnsatzSpace/FunctionEvaluator.hpp"
#include <iostream>

int main() {
  int refinement_level = 4;
  int polynomial_degree = 0;
  Bembel::Geometry geometry("sphere.dat");
  std::cout << "The geometry has " << geometry.get_number_of_patches()
            << " patches." << std::endl;
  Bembel::AnsatzSpace<Bembel::DummyOperator>(geometry, refinement_level,
                                             polynomial_degree);


  std::cout << Bembel::LinearOperatorTraits<Bembel::DummyOperator>::OperatorOrder << std::endl;
  return 0;
}
