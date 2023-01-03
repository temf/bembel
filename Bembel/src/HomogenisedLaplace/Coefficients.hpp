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

#ifndef PI
#define PI 3.141592653589793
#endif

#include <functional>

#include <Eigen/Dense>
#include <Bembel/Quadrature>

namespace Bembel {

using namespace Eigen;

VectorXd getCoefficients(double precision);

inline unsigned int getDegree(double precision);

VectorXd getDisplacement(MatrixXd ps_l, MatrixXd ps_f, MatrixXd ps_b);

inline double k_mod(Vector3d in);

inline Vector3d Dk_mod(Vector3d in);


VectorXd getCoefficients(double precision) {

	VectorXd coeffs;

	unsigned int deg = getDegree(precision);
	unsigned int Msquare = POINT_DEGREE*POINT_DEGREE;

	 GaussSquare<POINT_DEGREE> GS;
	 MatrixXd xs = GS[POINT_DEGREE].xi_ - 0.5;

	 Vector3d ex(1.0, 0.0, 0.0);
	 Vector3d ey(0.0, 1.0, 0.0);
	 Vector3d ez(0.0, 0.0, 1.0);

	 MatrixXd ps_left(3, xs.cols());
	 ps_left.row(0) = -0.5*Eigen::VectorXd::Ones(xs.cols());
	 ps_left.block(1, 0, 2, xs.cols()) = xs.block(0, 0, 2, xs.cols());

	 MatrixXd ps_front(3, xs.cols());
	 ps_front.row(0) = xs.row(0);
	 ps_front.row(1) = -0.5*Eigen::VectorXd::Ones(xs.cols());
	 ps_front.row(2) = xs.row(1);

	 MatrixXd ps_bottom(3, xs.cols());
	 ps_bottom.block(0, 0, 2, xs.cols()) = xs.block(0, 0, 2, xs.cols());
	 ps_bottom.row(2) = -0.5*Eigen::VectorXd::Ones(xs.cols());

	 VectorXd displacement = getDisplacement(ps_left, ps_front, ps_bottom);

	 MatrixXd systemMatrix(6*Msquare, ((deg+1)*(deg+2))/2);

	 return coeffs;


}


/**
 * Returns the degree n of the sphericals expansion
 * given a precision. Compare the Plots in the msc thesis.
 */
inline unsigned int getDegree(double precision) {
	if(precision > 1e-4) {
		return 4;
	} else if(precision > 1e-6) {
		return 8;
	} else {
		return 12;
	}
}

/**
 * Returns the displacement field for the homogenised
 * Laplace calculation. This corresponds to the
 * right-hand side of the system
 */
VectorXd getDisplacement(MatrixXd ps_l, MatrixXd ps_f, MatrixXd ps_b) {

	std::function<double(Vector3d)>  u = [](Vector3d in) {return k_mod(in);};
	std::function<Vector3d(Vector3d)> Du = [](Vector3d in) {return Dk_mod(in);};

	unsigned int Msquare = POINT_DEGREE*POINT_DEGREE;
	unsigned int k;

	Vector3d ex(1.0, 0.0, 0.0);
	Vector3d ey(0.0, 1.0, 0.0);
	Vector3d ez(0.0, 0.0, 1.0);

	Vector3d tmp;
	VectorXd d(6*Msquare); /* the return vector */

	/* left - right displacement */
	for(k = 0; k < Msquare; k++) {
		tmp = Du(ps_l.col(k)) - Du(ps_l.col(k) + ex);
		d(k) = u(ps_l.col(k)) - u(ps_l.col(k) + ex);
		d(k + Msquare) = tmp(0);
	}

	/* front - back displacement */
	for(k = 0; k < Msquare; k++) {
		tmp = Du(ps_f.col(k)) - Du(ps_f.col(k) + ey);
		d(k + 2*Msquare) = u(ps_f.col(k)) - u(ps_f.col(k) + ey);
		d(k + 3*Msquare) = tmp(1);
	}

	/* bottom - top displacement */
	for(k = 0; k < Msquare; k++) {
		tmp = Du(ps_b.col(k)) - Du(ps_b.col(k) + ez);
		d(k + 4*Msquare) = u(ps_b.col(k)) - u(ps_b.col(k) + ez);
		d(k + 5*Msquare) = tmp(2);
	}

	return d;

}

inline double k_mod(Vector3d in) {

	double r = 0.0;

	short i, j, k;
	Vector3d m;

	for(i = -1; i <= 1; i++) {
		for(j = -1; j <= 1; j++) {
			for(k = -1; k <= 1; k++) {
				m = Vector3d(i, j, k);
				r += 1.0/((in - m).norm());
			}
		}
	}

	r /= (4*PI);

	/* the part to ensure the vanishing mean on the Laplacian */
	r += (in.dot(in))/6.0;

	return r;
}

inline Vector3d Dk_mod(Vector3d in) {

	Vector3d r, s;
	double snorm;
	r.setZero();

	short i, j, k;
	for(i = -1; i <= 1; i++) {
		for(j = -1; j <= 1; j++) {
			for(k = -1; k <= 1; k++) {
				s = in - Vector3d(i, j, k);
				snorm = s.norm();
				r -= s/(snorm*snorm*snorm);
			}
		}
	}

	r /=(4*PI);

	/* the part to ensure the vanishing mean on the Laplacian */
	r += in/3.0;


	return r;
}

} /* namespace Bembel */





#endif /* BEMBEL_HOMOGENISEDLAPLACE_COEFFICIENTS_H_ */
