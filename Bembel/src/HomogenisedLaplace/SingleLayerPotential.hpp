/*
 * SingleLayerPotential.hpp
 *
 *  Created on: 9 Jan 2023
 *      Author: ricrem00
 */

#ifndef BEMBEL_LINEAROPERATOR_HOMOGENISEDLAPLACE_HOMOGENISEDLAPLACESINGLELAYERPOTENTIAL_H_
#define BEMBEL_LINEAROPERATOR_HOMOGENISEDLAPLACE_HOMOGENISEDLAPLACESINGLELAYERPOTENTIAL_H_

namespace Bembel {
// forward declaration of class HomogenisedLaplaceSingleLayerPotential in order to define
// traits
template <typename LinOp>
class HomogenisedLaplaceSingleLayerPotential;

template <typename LinOp>
struct PotentialTraits<HomogenisedLaplaceSingleLayerPotential<LinOp>> {
	typedef Eigen::VectorXd::Scalar Scalar;
	static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup HomogenisedLaplace
 */
template <typename LinOp>
class HomogenisedLaplaceSingleLayerPotential
		: public PotentialBase<HomogenisedLaplaceSingleLayerPotential<LinOp>, LinOp> {
	// implementation of the kernel evaluation, which may be based on the
	// information available from the superSpace

		private:
	unsigned int deg;
	Eigen::VectorXd cs;

		public:
	HomogenisedLaplaceSingleLayerPotential() {
		this->deg = getDegree(HomogenisedLaplaceSingleLayerOperator::getPrecision());
		this->cs = getCoefficients(HomogenisedLaplaceSingleLayerOperator::getPrecision());
	}
	Eigen::Matrix<
	typename PotentialReturnScalar<
	typename LinearOperatorTraits<LinOp>::Scalar, double>::Scalar,
	1, 1>
	evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
			const ElementTreeNode &element,
			const Eigen::Vector3d &point,
			const SurfacePoint &p) const {
		// get evaluation points on unit square
		auto s = p.segment<2>(0);

		// get quadrature weights
		auto ws = p(2);

		// get points on geometry and tangential derivatives
		auto x_f = p.segment<3>(3);
		auto x_f_dx = p.segment<3>(6);
		auto x_f_dy = p.segment<3>(9);

		// compute surface measures from tangential derivatives
		auto x_kappa = x_f_dx.cross(x_f_dy).norm();

		// evaluate kernel
		auto kernel = evaluateKernel(point, x_f);

		// assemble Galerkin solution
		auto cauchy_value = fun_ev.evaluate(element, p);

		// integrand without basis functions
		auto integrand = kernel * cauchy_value * x_kappa * ws;

		return integrand;
	}

	/**
	 * \brief Fundamental solution of the homogenised Laplace problem
	 */
	double evaluateKernel(const Eigen::Vector3d &x,
			const Eigen::Vector3d &y) const {

		return k_mod(x-y) + evaluate_solid_sphericals(x-y, this->cs, this->deg, false);
	}

};

}  // namespace Bembel




#endif /* BEMBEL_LINEAROPERATOR_HOMOGENISEDLAPLACE_HOMOGENISEDLAPLACESINGLELAYERPOTENTIAL_H_ */
