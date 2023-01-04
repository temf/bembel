/*
 * Sphericals.hpp
 *
 *  Created on: 3 Jan 2023
 *      Author: remo
 */

#ifndef BEMBEL_UTIL_SPHERICALS_H_
#define BEMBEL_UTIL_SPHERICALS_H_
#include <Eigen/Dense>

#if !defined pi
#define pi 3.1415926535897932385
#endif

namespace Bembel {

double evaluate_sphericals(Eigen::Vector3d x, Eigen::VectorXd k, unsigned int nk);

Eigen::Vector3d evaluate_dsphericals(Eigen::Vector3d x, Eigen::VectorXd k, unsigned int nk);

Eigen::Matrix<double, Eigen::Dynamic, 2> spherical_harmonics_full(Eigen::Vector3d x, unsigned int N);

inline Eigen::Vector2d spherical_prev(Eigen::Vector3d x, int m, int n, Eigen::Vector2d y1, Eigen::Vector2d y2);

Eigen::Matrix<double, 3, Eigen::Dynamic> Dsolid_harmonics_full(Eigen::Vector3d x, unsigned int N, Eigen::VectorXd L);

Eigen::Matrix<double, 3, 2> dsolid_spherical_prev(Eigen::Vector3d y, unsigned int m, unsigned int n, Eigen::VectorXd L, double y_re, double y_im);

Eigen::VectorXd legendreFull(unsigned int N, double t);

double constant(unsigned int m, unsigned int n);

Eigen::Matrix3d functionalMatrix(Eigen::Vector3d z);

/**
 * Evaluates the series \sum_{n = 0}^{nk} \sum_{m = -n}^n k_m^n Y_m^n(x) for real coefficients k
 *
 * Input:	x:		Eigen::Vector3d		The point of evaluation, a vector with length 1
 * 			k:		Eigen::VectorXd		The coefficients stored in the order [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)]
 * 			nk:		unsigned int		The degree
 */
double evaluate_sphericals(Eigen::Vector3d x, Eigen::VectorXd k, unsigned int nk) {
	unsigned int	m, n;
	double		z1[2], z2[2], z3[2], z1_start[2];
	double		r, fac, rootTimesZ, root_2, root_3;

	r = z1[1] = 0;
	z1[0] = 0.5/sqrt(pi);
	if(nk <= 1)	{return k(0)*z1[0];}

	for (m=0; m<nk-1; m++) {
		if(m == 0)	{fac = 1.0;}
		else		{fac = 2.0;}

		z1_start[0] = z1[0];
		z1_start[1] = z1[1];
		rootTimesZ = sqrt(2*m+3)*x(2);
		z2[0] = rootTimesZ*z1[0];
		z2[1] = rootTimesZ*z1[1];
		r += fac*k( m   *(m+1)+m)*z1[0]; // + k[ m   *(m+1)-m]*z1[1];
		r += fac*k((m+1)*(m+2)+m)*z2[0]; // + k[(m+1)*(m+2)-m]*z2[1];
		for (n=m+2; n<nk; n++) {
			root_2 = sqrt((2*n+1.0)/((n-m)*(n+m)));
			rootTimesZ = sqrt(2*n-1)*x(2);
			root_3 = sqrt(((n+m-1.0)*(n-m-1.0))/(2*n-3));
			z3[0] = root_2*(rootTimesZ*z2[0] - root_3*z1[0]);
			z3[1] = root_2*(rootTimesZ*z2[1] - root_3*z1[1]);
			r += fac*k(n*(n+1)+m)*z3[0]; // + k[n*(n+1)-m]*z3[1];
			z1[0] = z2[0];
			z1[1] = z2[1];
			z2[0] = z3[0];
			z2[1] = z3[1];
		}
		root_2 = sqrt((2*m+3.0)/(2*m+2));
		z1[0] = root_2*(x(0)*z1_start[0]-x(1)*z1_start[1]);
		z1[1] = root_2*(x(0)*z1_start[1]+x(1)*z1_start[0]);
	}
	r += 2*k((nk-1)*(nk+1))*z1[0]; //+k[(nk-1)*(nk-1)]*z1[1];
	return r;
}

/**
 * Evaluates the series \sum_{n = 0}^{nk} \sum_{m = -n}^n k_m^n  grad Y_m^n(x) for real coefficients k
 *
 * Input:	x:		Eigen::Vector3d		The point of evaluation, a vector with length 1
 * 			k:		Eigen::VectorXd		The coefficients stored in the order [(0,  0), (1, -1), (1, 0), (1, 1)
 * 																			  (2, -2), (2, -1), ... (n, n)]
 * 			nk:		unsigned int		The degree
 */
Eigen::Vector3d evaluate_dsphericals(Eigen::Vector3d x, Eigen::VectorXd k, unsigned int nk) {
	unsigned int		m, n;
	double				z1[2], z2[2], z3[2], z1_start[2];
	Eigen::Vector3d		dr;
	double 				fac, rootTimesZ, root_2, root_3;

	z1[0] = sqrt(0.375/pi);
	z1[1] = 0;

	dr = Eigen::Vector3d(0.0, 0.0, 0.0);

	if(nk == 1)	{return dr;}

	for (m=1; m<nk-1; m++) {
		z1_start[0] = z1[0];
		z1_start[1] = z1[1];
		rootTimesZ = sqrt(2*m+3)*x(2);
		z2[0] = rootTimesZ*z1[0];
		z2[1] = rootTimesZ*z1[1];
		dr(0) += 2*m*(k( m   *(m+1)+m)*z1[0]); //+a[ m   *(m+1)-m]*z1[1]);
		dr(0) += 2*m*(k((m+1)*(m+2)+m)*z2[0]); //+a[(m+1)*(m+2)-m]*z2[1]);
		dr(1) -= 2*m*(k( m   *(m+1)+m)*z1[1]); //-a[ m   *(m+1)-m]*z1[0]);
		dr(1) -= 2*m*(k((m+1)*(m+2)+m)*z2[1]); //-a[(m+1)*(m+2)-m]*z2[0]);

		if(m == 1) {
			fac = 1.0;
			dr(2) += fac*sqrt(2*m  )*(k( m   *(m+1)+(m-1))*z1[0]+k( m   *(m+1)-(m-1))*z1[1]);
			dr(2) += fac*sqrt(4*m+2)*(k((m+1)*(m+2)+(m-1))*z2[0]+k((m+1)*(m+2)-(m-1))*z2[1]);
		} else {
			fac = 2.0;
			dr(2) += fac*sqrt(2*m  )*(k( m   *(m+1)+(m-1))*z1[0]); //+a[ m   *(m+1)-(m-1)]*z1[1]);
			dr(2) += fac*sqrt(4*m+2)*(k((m+1)*(m+2)+(m-1))*z2[0]); //+a[(m+1)*(m+2)-(m-1)]*z2[1]);
		}

		for (n=m+2; n<nk; n++) {
			root_2 = sqrt((2*n+1.0)/((n-m)*(n+m)));
			rootTimesZ = sqrt(2*n-1)*x(2);
			root_3 = sqrt(((n+m-1.0)*(n-m-1.0))/(2*n-3));
			z3[0] = root_2*(rootTimesZ*z2[0] - root_3*z1[0]);
			z3[1] = root_2*(rootTimesZ*z2[1] - root_3*z1[1]);
			dr(0) += 2*m*(k(n*(n+1)+m)*z3[0]); //+a[n*(n+1)-m]*z3[1]);
			dr(1) -= 2*m*(k(n*(n+1)+m)*z3[1]); //-a[n*(n+1)-m]*z3[0]);
			dr(2) += fac*sqrt((n+m)*(n-m+1))*(k(n*(n+1)+(m-1))*z3[0]); //+a[n*(n+1)-(m-1)]*z3[1]);
			z1[0] = z2[0];
			z1[1] = z2[1];
			z2[0] = z3[0];
			z2[1] = z3[1];
		}
		root_2 = sqrt((2*m+3.0)/(2*m+2));
		z1[0] = root_2*(x(0)*z1_start[0]-x(1)*z1_start[1]);
		z1[1] = root_2*(x(0)*z1_start[1]+x(1)*z1_start[0]);
	}
	dr(0) += (nk-1)*(2*k((nk-1)*(nk+1))*z1[0]); //+a[(na-1)*(na-1)]*z1[1]);
	dr(1) -= (nk-1)*(2*k((nk-1)*(nk+1))*z1[1]); //-a[(na-1)*(na-1)]*z1[0]);
	dr(2) += 2*sqrt(2*(nk-1))*(k(nk*(nk-1)+(nk-2))*z1[0]+k(nk*(nk-1)-(nk-2))*z1[1]);

	return dr;
}

/**
 * Calculates the the solid harmonics phi_n^m(x), ordered by [phi_0^0, phi_1^{-1}, phi_1^0 phi_1^1, phi_2^{-2}, ..., phi_N^N]
 */
Eigen::Matrix<double, Eigen::Dynamic, 2> spherical_harmonics_full(Eigen::Vector3d x, unsigned int N) {

	Eigen::VectorXd real((N+1)*(N+1));
	Eigen::VectorXd imag((N+1)*(N+1));

	Eigen::Vector3d y = x/x.norm();

	int m, n;

	/* handle n = 0 separately */
	real(0) = 0.5/sqrt(pi);
	imag(0) = 0.0;

	/* assemble two temporary vectors */
	Eigen::Vector2d z1, z2, tmp;

	for(n=1; n <= N; n++) {

		for(m=0; m < n-1; m++) {
			z1(0) = real((n-1)*(n-1) + n-1 + m);
			z1(1) = real((n-2)*(n-2) + n-2 + m);
			z2(0) = imag((n-1)*(n-1) + n-1 + m);
			z2(1) = imag((n-2)*(n-2) + n-2 + m);

			tmp = spherical_prev(y, m, n, z1, z2);

			real(n*n + n+m) = tmp(0);
			imag(n*n + n+m) = tmp(1);

			// if m > 0, copy the value with the prefactor
			if(m > 0) {
				real(n*n + n-m) = tmp(0);
				imag(n*n + n-m) = -tmp(1);
			}
		}


		/* for m = n-1 */
		m = n-1;
		z1(0) = real(m*m + m+m);
		z2(0) = imag(m*m + m+m);

		tmp = spherical_prev(y, m, n, z1, z2);
		real(n*n + n+m) = tmp(0);
		imag(n*n + n+m) = tmp(1);

		if(m > 0) {
			real(n*n + n-m) = tmp(0);
			imag(n*n + n-m) = -tmp(1);
		}


		/* for m = n */
		m = n;
		z1(0) = real((n-1)*(n-1) + n-1 + n-1);
		z2(0) = imag((n-1)*(n-1) + n-1 + n-1);

		tmp = spherical_prev(y, m, n, z1, z2);
		real(n*n + n+m) = tmp(0);
		imag(n*n + n+m) = tmp(1);

		if(m > 0) {
			real(n*n + n-m) = tmp(0);
			imag(n*n + n-m) = -tmp(1);
		}


	}

	Eigen::Matrix<double, (N+1)*(N+1), 2> res;
	res.col(0) = real;
	res.col(1) = imag;

	return res;


}

/*
 * Calculates the spherical harmonics value Y_n^m(x) based on the previous values
 */
inline Eigen::Vector2d spherical_prev(Eigen::Vector3d x, int m, int n, Eigen::Vector2d y1, Eigen::Vector2d y2) {
	Eigen::Vector2d z;

	if ((m == 0) && (n == 0)) {
		z(0) = 0.5/sqrt(pi);
		z(1) = 0;
	} else if (m == n) {
		z(0) = sqrt((2*m+1.0)/(2*m))*(x(0)*y1(0)-x(1)*y2(0));
		z(1) = sqrt((2*m+1.0)/(2*m))*(x(0)*y2(0)+x(1)*y1(0));
	} else if (m+1 == n) {
		z(0) = y1(0)*sqrt(2*m+3)*x(2);
		z(1) = y2(0)*sqrt(2*m+3)*x(2);
	} else {
		z(0) = sqrt((2*n+1.0)/((n-m)*(n+m)))*(sqrt(2*n-1)*x(2)*y1(0) - sqrt(((n+m-1.0)*(n-m-1.0))/(2*n-3))*y1(1));
		z(1) = sqrt((2*n+1.0)/((n-m)*(n+m)))*(sqrt(2*n-1)*x(2)*y2(0) - sqrt(((n+m-1.0)*(n-m-1.0))/(2*n-3))*y2(1));
	}

	return z;
}

/**
 * Calculates all gradients of the solid harmonics, given the Legendre Coefficients L
 */
Eigen::Matrix<double, 3, Eigen::Dynamic> Dsolid_harmonics_full(Eigen::Vector3d x, unsigned int N, Eigen::MatrixXd ys) {

	Eigen::Matrix<double, 3, Eigen::Dynamic> reals, imags;

	Eigen::Vector3d y = x/x.norm();

	Eigen::VectorXd L = legendreFull(N, y(2));

	Eigen::Matrix<double, 3, 2> z;

	int m, n;

	for(n = 0; n <= N; n++) {
		for(m = 0; m <= n; m++) {
			z = dsolid_spherical_prev(y, m, n, L, ys(n*n + n+m, 0), ys(n*n + n+m, 1));

			reals(n*n + n+m, 0) = z(0, 0);
			reals(n*n + n+m, 1) = z(1, 0);
			reals(n*n + n+m, 2) = z(2, 0);

			imags(n*n + n+m, 0) = z(0, 1);
			imags(n*n + n+m, 1) = z(1, 1);
			imags(n*n + n+m, 2) = z(2, 1);


			if(m > 0) {
				reals(n*n + n-m, 0) = z(0, 0);
				reals(n*n + n-m, 1) = z(1, 0);
				reals(n*n + n-m, 2) = z(2, 0);

				imags(n*n + n-m, 0) = -z(0, 1);
				imags(n*n + n-m, 1) = -z(1, 1);
				imags(n*n + n-m, 2) = -z(2, 1);
			}

		}
	}

	return reals;

}

inline Eigen::Matrix<double, 3, 2> dspherical_prev(Eigen::Vector3d x, unsigned int m, unsigned int n, Eigen::VectorXd L) {

	double		c;
	unsigned int	i;

	Eigen::Matrix<double, 3, 2> z;

	if (m == 0) {
		z(0, 0) = z(0, 1) = z(2, 1) = 0;
		if(1 > n) {
			z(2, 0) = 0;
		} else {
			z(2, 0) = L((n*(n+1))/2 + 1);
		}
	} else {
		if(m > n) {
			z(0, 0) = 0;
		} else {
			z(0, 0) = L((n*(n+1))/2 + m);
		}

		if(m+1 > n) {
			z(2, 0) = 0;
		} else {
			z(2, 0) = L((n*(n+1))/2 + m+1);
		}
		z(0, 1) = z(2, 1) = 0;

		c		= x(0) * z(2, 0) - x(1) * z(2, 1);
		z(2, 1) = x(0) * z(2, 1) + x(1) * z(2, 0);
		z(2, 0) = c;
		for (i=1; i<m; i++) {
			c     = x(0) * z(0, 0) - x(1) * z(0, 1);
			z(0, 1) = x(0) * z(0, 1) + x(1) * z(0, 0);
			z(0, 0) = c;

			c     = x(0) * z(2, 0) - x(1) * z(2, 1);
			z(2, 1) = x(0) * z(2, 1) + x(1) * z(2, 0);
			z(2, 0) = c;
		}
	}

	c = constant(m,n);
	z(0, 0) *= m*c;
	z(0, 1) *= m*c;
	z(1, 0) = -z(0, 1);
	z(1, 1) = +z(0, 0);
	z(2, 0) *= c;
	z(2, 1) *= c;

	return z;
}

/**
 * Calculates the derivative of the solid harmonics function phi_n^m at the point y with the spherical values
 */
Eigen::Matrix<double, 3, 2> dsolid_spherical_prev(Eigen::Vector3d y, unsigned int m, unsigned int n, Eigen::VectorXd L, double y_re, double y_im) {

	//normalise y
	// x denotes the normalised y and is passed to the original functions

	double r = y.norm();

	Eigen::Vector3d x = y/r;

	Eigen::Matrix<double, 3, 2> z_harm;

	Eigen::Matrix3d A;

	unsigned int m_tilde;

	// part with the gradient of Y
	if(m >= 0) {

		m_tilde = m;

		z_harm = dspherical_prev(x, m_tilde, n, L);

		// multiply the gradient of Y by r^(n-3) and the matrix A
		// pow has a problem if the exponent is negative, so we adapt it here
		double r_n3;
		if(n > 3) {
			r_n3 = pow(r, n-3);
		} else if (n == 3) {
			r_n3 = 1.0;
		} else {
			r_n3 = pow(1/r, 3-n);
		}

		A = functionalMatrix(y);
		z_harm = r_n3*(A*z_harm);

	} else {
		m_tilde = -m;
		// define r^(n-3) as above
		double r_n3;
		if(n > 3) {
			r_n3 = pow(r, n-3);
		} else if (n == 3) {
			r_n3 = 1.0;
		} else {
			r_n3 = pow(1/r, 3-n);
		}

		// get the part of the gradient of Y, conjugate and multiply by the functional matrix
		z_harm = dspherical_prev(x, m_tilde, n, L);
		A = functionalMatrix(y);
		z_harm = r_n3*(A*z_harm);

		// get the -1^m and complex conjugation

		z_harm.col(1) *= -1;

	}

	// the radial part just needs an adaptation for the Y of the spherical function
	Eigen::Matrix<double, 3, 2> z_rad;

	// case n = 0: r^n = 1 is constant
	if(n == 0) {
		z_rad.setZero();
	} else {
		// the gradient of r^n is n*r^(n-2)*y = nr^(n-1)*x
		// again adapt the power
		double nr_n1;
		if(n > 1) {
			nr_n1= n*pow(r, n-1);
		} else if(n == 1) {
			nr_n1 = 1;
		} else {
			nr_n1 = n*pow(1/r, 1-n);
		}


		double y_imt;
		if(m >= 0)	{y_imt = y_im;}
		else		{y_imt = -y_im;}

		z_rad.col(0) =  y_re*nr_n1*x;
		z_rad.col(1) = y_imt*nr_n1*x;
	}

	return z_harm + z_rad;

}


/**
 * Returns the values of the spherical polynomials P_n^m(t)
 */
Eigen::VectorXd legendreFull(unsigned int N, double t) {

	Eigen::VectorXd L(((N+1)*(N+2))/2);

	int n, s, m;

	/* n = 0 */
	L(0) = 1;

	/* n = 1 */
	L(1) = t;
	L(2) = 1;

	for(n = 2; n <= N; n++) {
		s = (n*(n+1))/2;
		for(m = 0; m < n-1; m++) {
			L(s+m) = ((2*n-1)*t*L((n*(n-1))/2 + m) - (n+m-1)*L(((n-1)*(n-2))/2 + m))/(n-m);
		}

		//m = n-1
		m = n-1;
		L(s+m) = (2*m + 1)*t*L((m*(m+1))/2 + m);

		//m = n
		m = n;
		L(s+m) = (2*m - 1)*L((n*(n-1))/2 + n-1);
	}

	return L;
}

/**
 *  gives the norming factor of the spherical polynomials
 */
double constant(unsigned int m, unsigned int n) {
	unsigned int 	i;
	double		c;
	c = sqrt((2*n+1)/(4*pi));
	for (i=0; i<m; i++)	{c *= 1/sqrt((n+i+1)*(n-i));}

	return c;
}

/**
 * Gives the Jacobi Matrix of the transformation z |--> z / |z|
 */
Eigen::Matrix3d functionalMatrix(Eigen::Vector3d z) {

		Eigen::Matrix3d M;

		double r = z.norm();

		M(0, 0) = r*r - z(0)*z(0);
		M(1, 1) = r*r - z(1)*z(1);
		M(2, 2) = r*r - z(2)*z(2);


		M(1, 0) = M(0, 1) = -z(0)*z(1);
		M(2, 0) = M(0, 2) = -z(0)*z(2);
		M(2, 1) = M(1, 2) = -z(1)*z(2);

		return M;

}


} // namespace Bembel


#endif /* BEMBEL_UTIL_SPHERICALS_H_ */
