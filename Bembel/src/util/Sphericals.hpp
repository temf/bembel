/*
 * Sphericals.hpp
 *
 *  Created on: 18.03.2021
 *      Author: remo
 */

#ifndef BEMBEL_UTIL_SPHERICALS_H_
#define BEMBEL_UTIL_SPHERICALS_H_
#include <Eigen/Dense>

#if !defined pi
#define pi 3.1415926535897932385
#endif

namespace Bembel {

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



} // namespace Bembel


#endif /* BEMBEL_UTIL_SPHERICALS_H_ */
