/*
 * test_HomogenisedCoefficients.cpp
 *
 *  Created on: 5 Jan 2023
 *      Author: ricrem00
 */


#include <iostream>
#include <functional>

#include <math.h>

#include <Bembel/HomogenisedLaplace>
#include <Eigen/Dense>

inline double k_per(Eigen:: Vector3d in, Eigen::VectorXd coeffs, unsigned int deg);

inline Eigen::Vector3d Dk_per(Eigen:: Vector3d in, Eigen::VectorXd coeffs, unsigned int deg);

Eigen::Matrix<double, 2, Eigen::Dynamic> tensorise(Eigen::ArrayXd xs);

int main() {

	using namespace Bembel;
	using namespace Eigen;

	double precision = 1e-12;
	unsigned int deg = getDegree(precision);

	VectorXd cs = getCoefficients(precision);

	unsigned int Npoints = 101;

	ArrayXd xs = ArrayXd::LinSpaced(Npoints, -0.5, 0.5);
	MatrixXd ys = tensorise(xs);

	double err = 0.0;

	Vector3d ex(1.0, 0.0, 0.0);
	Vector3d ey(0.0, 1.0, 0.0);
	Vector3d ez(0.0, 0.0, 1.0);

	std::function<double(Vector3d)> u = [cs, deg](Vector3d in) 									{return k_per(in, cs, deg);};
	std::function<double(Vector3d, unsigned int)> Du = [cs, deg](Vector3d in, unsigned int d)	{return Dk_per(in, cs, deg)(d);};


	unsigned int k;
	Vector3d v;
	for(k = 0; k < Npoints*Npoints; k++) {
		v = Vector3d(-0.5, ys(0, k), ys(1, k));
		err += fabs(k_per(v, cs, deg) - k_per(v + ex, cs, deg));
		err += fabs(Du(v, 0) - Du(v + ex, 0));

		v = Vector3d(ys(0, k), -0.5, ys(1, k));
		err += fabs(u(v) - u(v + ey));
		err += fabs(Du(v, 1) - Du(v + ey, 1));

		v = Vector3d(ys(0, k), ys(1, k), -0.5);
		err += fabs(u(v) - u(v + ez));
		err += fabs(Du(v, 2) - Du(v + ez, 2));
	}

	err /= Npoints*Npoints;

	std::cout << "Average Pointwise Error is " << err << std::endl;

	return 0;
}

inline double k_per(Eigen::Vector3d in, Eigen::VectorXd coeffs, unsigned int deg) {

	double norm = in.norm();
	double rn;


	Eigen::VectorXd coeffs_p(coeffs.rows());
	unsigned int n;

	/* add the radial weights */
	coeffs_p(0) = coeffs(0);
	coeffs_p.segment(1, 3) = norm*coeffs.segment(1, 3);
	rn = norm;
	for(n = 2; n <= deg; n++) {
		rn *= norm;
		coeffs_p.segment(n*n, 2*n+1) = rn*coeffs.segment(n*n, 2*n+1);
	}

	return Bembel::k_mod(in) + Bembel::evaluate_sphericals(in/norm, coeffs_p, deg);
}

inline Eigen::Vector3d Dk_per(Eigen::Vector3d in, Eigen::VectorXd coeffs, unsigned int deg) {

	if(deg == 0) {return Eigen::Vector3d(0.0, 0.0, 0.0);}

	double norm = in.norm();
	double nr_n1, r_n3;
	unsigned int n, m;

	Eigen::Vector3d y = in/norm;
	Eigen::VectorXd c_harm(coeffs.rows()), c_rad(coeffs.rows());

	for(n = 0; n <= deg; n++) {
		if(n == 0)	{
			c_harm(0) = 0;
			c_rad(0) = 0;
			continue;
		} else if(n == 1) {
			nr_n1 = 1.0;
			r_n3 = 1.0/(norm*norm);
		} else if(n == 2) {
			nr_n1 = 2*norm;
			r_n3 = 1.0/norm;
		} else if(n == 3) {
			nr_n1 = 3*norm*norm;
			r_n3 = 1.0;
		} else {
			nr_n1 = n*pow(norm, n-1);
			r_n3 = pow(norm, n-3);
		}

		c_harm.segment(n*n, 2*n+1) = r_n3*coeffs.segment(n*n, 2*n+1);
		c_rad.segment(n*n, 2*n+1) = nr_n1*coeffs.segment(n*n, 2*n+1);
	}


	Eigen::Vector3d res = Bembel::Dk_mod(in);

	res += y*Bembel::evaluate_sphericals(y, c_rad, deg);
	res += Bembel::functionalMatrix(in)*Bembel::evaluate_dsphericals(y, c_harm, deg);

	return res;
}

Eigen::Matrix<double, 2, Eigen::Dynamic> tensorise(Eigen::ArrayXd xs) {

	unsigned int len = xs.rows();


	Eigen::MatrixXd ys(2, len*len);

	unsigned int k;
	for(k = 0; k < len; k++) {
		ys.block(0, k*len, 1, len) = xs(k)*Eigen::MatrixXd::Ones(1, len);
		ys.block(1, k*len, 1, len) = xs.transpose();
	}

	return ys;

}
