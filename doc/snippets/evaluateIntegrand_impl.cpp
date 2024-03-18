// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

// [operator]
void evaluateIntegrand_impl(
    const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
    Eigen::Matrix<
        typename LinearOperatorTraits<LaplaceSingleLayerOperator>::Scalar,
        Eigen::Dynamic, Eigen::Dynamic> *intval);
// [operator]

// [potential]
Eigen::Matrix<typename PotentialReturnScalar<
                  typename LinearOperatorTraits<LinOp>::Scalar, double>::Scalar,
              1, 1>
evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                       const ElementTreeNode &element,
                       const Eigen::Vector3d &point, const SurfacePoint &p);
// [potential]

// [linearform]
void evaluateIntegrand_impl(const T &super_space, const SurfacePoint &p,
                            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *intval);
// [linearform]
