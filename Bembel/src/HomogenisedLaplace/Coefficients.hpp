/*
 * Coefficients.hpp
 *
 *  Created on: 3 Jan 2023
 *      Author: ricrem00
 */

#ifndef BEMBEL_HOMOGENISEDLAPLACE_COEFFICIENTS_H_
#define BEMBEL_HOMOGENISEDLAPLACE_COEFFICIENTS_H_

#ifndef POINT_DEGREE
#define POINT_DEGREE 20 /* the number of points on the surface of the cube */
#endif

#include <Eigen/Dense>
#include <Bembel/Quadrature>

namespace Bembel{

Eigen::VectorXd getCoefficients(double precision);

unsigned int getDegree(double precision);


Eigen::VectorXd getCoefficients(double precision) {

	unsigned int deg = getDegree(precision);

	 GaussSquare<POINT_DEGREE> GS;
	 Eigen::MatrixXd xs = GS[POINT_DEGREE].xi_;


}


/**
 * Returns the degree n of the sphericals expansion
 * given a precision. Compare the Plots in the msc thesis.
 */
unsigned int getDegree(double precision) {
	if(precision > 1e-4) {
		return 4;
	} else if(precision > 1e-6) {
		return 8;
	} else {
		return 12;
	}
}

}





#endif /* BEMBEL_HOMOGENISEDLAPLACE_COEFFICIENTS_H_ */
