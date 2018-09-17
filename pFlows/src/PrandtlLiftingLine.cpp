#include "PrandtlLiftingLine.h"
/*////////////////////////////////////////////////////////////////////////////
PrandtlLiftingLine.cpp

An implementation of the classical Prandtl lifting line theory

Copyright 2017 HJA Bird

mFlow is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mFlow is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mFlow.  If not, see <http://www.gnu.org/licenses/>.
*/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include <Eigen/Dense>

#include <HBTK/Constants.h>
#include <HBTK/GaussLegendre.h>

double mFlow::PrandtlLiftingLine::lift_coefficient()
{
	assert(n_terms > 0);

	double k_1 = -1. / (HBTK::Constants::pi());
	double k_2 = 1. / (4 * HBTK::Constants::pi());
	auto c_1 = [&](int k, double theta) { 
		return k_1 * sin((2 * k + 1) * theta) / wing.chord(theta_to_y(theta));
	};
	auto c_2 = [&](int k, double theta) { 
		return k_2 * (2 * k + 1) * integral(k, theta); 
	};
	std::vector<double> collocation_points = generate_collocation_points();

	Eigen::MatrixXd mat;
	mat.resize(n_terms, n_terms);
	Eigen::VectorXd aoa, solution;
	aoa.resize(n_terms);

	for (int i = 0; i < n_terms; i++) {
		double y_i = collocation_points[i];
		for (int k = 0; k < n_terms; k++) {
			double theta = y_to_theta(y_i);
			mat(i, k) = c_1(k, theta) + c_2(k, theta);
		}
		aoa(i) = -incidence(y_i) + zero_lift_incidence(y_i);
	}
	solution = mat.lu().solve(aoa);
	
	auto kernal = [&](double y) {
		double theta = y_to_theta(y);
		double vort = 0;
		for (int k = 0; k < n_terms; k++) {
			vort += solution(k) * sin((2 * k + 1) * theta);
		}
		return vort;
	};
	auto quad = HBTK::gauss_legendre(25);
	quad.linear_remap(-wing.semispan(), wing.semispan());
	double cl = 2 * quad.integrate(kernal) / wing.area();
	return cl;
}

double mFlow::PrandtlLiftingLine::y_to_theta(double y) const
{
	double a = 2 * y / wing.span;
	return acos(a);
}

double mFlow::PrandtlLiftingLine::theta_to_y(double theta) 
{
	return wing.semispan() * cos(theta);
}

std::vector<double> mFlow::PrandtlLiftingLine::generate_collocation_points()
{
	std::vector<double> collocation_points(n_terms);
	const auto hpi = HBTK::Constants::pi() / 2;
	for (auto idx = 0; idx < n_terms; idx++) {
		collocation_points[idx] = 
			theta_to_y((hpi + idx * HBTK::Constants::pi()) / (2 * n_terms));
	}
	return collocation_points;
}

double mFlow::PrandtlLiftingLine::integral(int k, double theta) const
{
	// Using the singularity subtraction method and Glauert integral.
	const double c_denom = cos(theta);
	const double c_num = cos((2 * k + 1) * theta);
	auto integrand = [&](double theta_0) {
		return (cos((2 * k + 1) * theta_0) - c_num) / 
			(c_denom - cos(theta_0));
	};
	auto quad = HBTK::gauss_legendre(20);
	quad.linear_remap(0, HBTK::Constants::pi());
	return 2 * quad.integrate(integrand) / wing.span;
}
