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
#define PI M_PI
#endif

#include <functional>

#include <Eigen/Dense>
#include <Bembel/HomogenisedLaplace>
#include <Bembel/Quadrature>

namespace Bembel {

using namespace Eigen;

VectorXd getCoefficients(double precision);

inline unsigned int getDegree(double precision);

VectorXd getDisplacement(MatrixXd ps_l, MatrixXd ps_f, MatrixXd ps_b);

inline double k_mod(Vector3d in);

inline Vector3d Dk_mod(Vector3d in);


VectorXd getCoefficients(double precision) {

	VectorXd diff;
	RowVectorXd difft;
	Vector3d v;
	MatrixXd spherical_values_pre, spherical_values_pst, Dsolid_values_pre, Dsolid_values_pst;

	unsigned int deg = getDegree(precision);
	unsigned int Msquare = (POINT_DEGREE + 1)*(POINT_DEGREE + 1);

	unsigned int m, n, k;
	double scale, fac, norm;

	GaussSquare<POINT_DEGREE> GS;
	MatrixXd xs = GS[POINT_DEGREE].xi_;
	xs -= 0.5*MatrixXd::Ones(xs.rows(), xs.cols());

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

	/* TODO: Is column major, but here row-major would make more sense... */
	MatrixXd systemMatrix(6*Msquare, ((deg+1)*(deg+2))/2 -1);

	for(k = 0; k < Msquare; k++) {
		/* left - right difference */
		v = ps_left.col(k);
		norm = v.norm();

		spherical_values_pre = spherical_harmonics_full(v, 	 	deg);
		spherical_values_pst = spherical_harmonics_full(v + ex, deg);

		Dsolid_values_pre = Dsolid_harmonics_full(v, 	  deg, spherical_values_pre);
		Dsolid_values_pst = Dsolid_harmonics_full(v + ex, deg, spherical_values_pst);

		for(n = 1; n <= deg; n++) {
			scale = pow(norm, n);
			diff = scale*(spherical_values_pre.block(n*n + n, 0, n+1, 1) - spherical_values_pst.block(n*n + n, 0, n+1, 1));
			difft = Dsolid_values_pre.block(0, n*n + n, 1, n+1) - Dsolid_values_pst.block(0, n*n + n, 1, n+1);

			/* adapt the scaling */
			diff.segment(1, n) *= 2.0;
			difft.segment(1, n) *= 2.0;

			systemMatrix.block(k, 			(n*(n+1))/2 - 1, 1, n+1) = diff.transpose();
			systemMatrix.block(k + Msquare, (n*(n+1))/2 - 1, 1, n+1) = difft;
		}

		/* front - back difference */
		v = ps_front.col(k);
		norm = v.norm();

		spherical_values_pre = spherical_harmonics_full(v, 	 	deg);
		spherical_values_pst = spherical_harmonics_full(v + ey, deg);

		Dsolid_values_pre = Dsolid_harmonics_full(v, 	  deg, spherical_values_pre);
		Dsolid_values_pst = Dsolid_harmonics_full(v + ey, deg, spherical_values_pst);

		for(n = 1; n <= deg; n++) {
			scale = pow(norm, n);
			diff = scale*(spherical_values_pre.block(n*n + n, 0, n+1, 1) - spherical_values_pst.block(n*n + n, 0, n+1, 1));
			difft = Dsolid_values_pre.block(1, n*n + n, 1, n+1) - Dsolid_values_pst.block(1, n*n + n, 1, n+1);

			/* adapt the scaling */
			diff.segment(1, n) *= 2.0;
			difft.segment(1, n) *= 2.0;

			systemMatrix.block(k + 2*Msquare, (n*(n+1))/2 - 1, 1, n+1) = diff.transpose();
			systemMatrix.block(k + 3*Msquare, (n*(n+1))/2 - 1, 1, n+1) = difft;
		}

		/* bottom - top difference */
		v = ps_bottom.col(k);
		norm = v.norm();

		spherical_values_pre = spherical_harmonics_full(v, 	 	deg);
		spherical_values_pst = spherical_harmonics_full(v + ez, deg);

		Dsolid_values_pre = Dsolid_harmonics_full(v, 	  deg, spherical_values_pre);
		Dsolid_values_pst = Dsolid_harmonics_full(v + ez, deg, spherical_values_pst);

		for(n = 1; n <= deg; n++) {
			scale = pow(norm, n);
			diff = scale*(spherical_values_pre.block(n*n + n, 0, n+1, 1) - spherical_values_pst.block(n*n + n, 0, n+1, 1));
			difft = Dsolid_values_pre.block(2, n*n + n, 1, n+1) - Dsolid_values_pst.block(2, n*n + n, 1, n+1);

			/* adapt the scaling */
			diff.segment(1, n) *= 2.0;
			difft.segment(1, n) *= 2.0;

			systemMatrix.block(k + 4*Msquare, (n*(n+1))/2 - 1, 1, n+1) = diff.transpose();
			systemMatrix.block(k + 5*Msquare, (n*(n+1))/2 - 1, 1, n+1) = difft;
		}

	}

	/* solve the system */
	VectorXd coeffs(((deg+1)*(deg+2))/2);
	coeffs.segment(1, ((deg+1)*(deg+2))/2-1) = systemMatrix.colPivHouseholderQr().solve(-displacement);

	/* Copy the stuff into the full Coefficient list */
	VectorXd coeffs_full((deg+1)*(deg+1));
	coeffs_full(0) = 0;
	for(n = 1; n <= deg; n++) {
		coeffs_full(n*n + n) = coeffs((n*(n+1))/2);
		for(m = 1; m <= n; m++) {
			coeffs_full(n*n + n+m) = coeffs((n*(n+1))/2 + m);
			coeffs_full(n*n + n-m) = coeffs((n*(n+1))/2 + m);
		}
	}

	/* calculate the first coefficient */




	return coeffs_full;


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

	unsigned int Msquare = (POINT_DEGREE + 1)*(POINT_DEGREE + 1);
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

	short int i, j, k;
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

	short int i, j, k;
	for(i = -1; i <= 1; i++) {
		for(j = -1; j <= 1; j++) {
			for(k = -1; k <= 1; k++) {
				s = in - Vector3d(i, j, k);
				snorm = s.norm();
				r -= s/(snorm*snorm*snorm);
			}
		}
	}

	r /= (4.0*PI);

	/* the part to ensure the vanishing mean on the Laplacian */
	r += in/3.0;


	return r;
}

} /* namespace Bembel */





#endif /* BEMBEL_HOMOGENISEDLAPLACE_COEFFICIENTS_H_ */
