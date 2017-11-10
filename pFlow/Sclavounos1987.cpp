#include "stdafx.h"
#include "Sclavounos1987.h"
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.cpp

Equations from the paper "An unsteady lifting-line theory", P.D. Sclavounos,
1987, Journal of Engineering Mathmatics

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

#include <array>

#include "../../HBTK/HBTK/Constants.h"
#include "../../HBTK/HBTK/GaussLegendre.h"
#include "../../HBTK/HBTK/Remaps.h"
#include "../../HBTK/HBTK/Integrators.h"

#include <Eigen/LU>

#include "Common.h"

namespace mFlow {
	Sclavounos1987::Sclavounos1987()
	{
	}


	Sclavounos1987::~Sclavounos1987()
	{
	}

	void Sclavounos1987::compute_collocation_points()
	{
		// See just after Eq7.1
		m_collocation_points.resize(number_of_terms);
		const double qpi = Constants::pi() / 4;
		for (int idx = 0; idx < number_of_terms; idx++) {
			m_collocation_points[idx] = qpi + idx * Constants::pi();
		}
		return;
	}

	std::complex<double> Sclavounos1987::d_3(double y)
	{
		double l; // semichord
		double v; // Normalised frequency
		std::complex<double> numerator, denominator;

		l = wing.chord_length(y);
		v = omega / U;

		numerator = 4 * U * exp(-Constants::i() * l * v);
		denominator = Constants::i() * Common::Hankle2_0(v * l) + Common::Hankle2_1(v * l);
		return numerator / denominator;
	}

	std::complex<double> Sclavounos1987::d_5(double y)
	{
		return -0.5 * wing.chord_length(y) * d_3(y);
	}

	std::complex<double> Sclavounos1987::K(double y)
	{
		std::complex<double> term_11, term_12, 
			term_121, term_122, term_123;
		double v;
		v = omega / U;

		term_11 = 0.5 * (y >= 0 ? 1. : -1.);

		term_121 = exp(-v * abs(y)) / abs(y);
		term_122 = Constants::i() * v * Ei(v * abs(y));
		term_123 = v * P(v * abs(y));
		term_12 = term_121 - term_122 + term_123;
		
		return term_11 * term_12;
	}

	void Sclavounos1987::compute_solution()
	{
		// Checks.
		assert(U > 0);
		assert(omega >= 0);
		assert(j == 3 || j == 5);
		assert(number_of_terms > 0);

		compute_collocation_points();

		// Allocation of matrices.
		m_solution.resize(number_of_terms);
		Eigen::MatrixXd gamma_matrix; // The matrix that directly gives vorticities at given points.
		// The matrix that representsthe integro-differential part.
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> integ_diff_matrix;
		// LHS matrix - the one we use for the solution (Yah, unneeded extra copies. I don't care).
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> LHS_matrix;
		// RHS vector
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> RHS_vector;
		gamma_matrix.resize(number_of_terms, number_of_terms);
		integ_diff_matrix.resize(number_of_terms, number_of_terms);
		LHS_matrix.resize(number_of_terms, number_of_terms);
		RHS_vector.resize(number_of_terms);

		// Compute gamma_matrix;
		for (int i = 0; i < number_of_terms; i++) {
			for (int j = 0; j < number_of_terms; j++) {
				gamma_matrix(i, j) = sin((2 * j) * m_collocation_points[i]);
			}
		}

		// The derivative of a fourier term with respect to y (ie: df/dtheta * dtheta/dy.
		auto sin_deriv_functor = [&](double y, int N)->double {
			// HJAB: see rough notes 09/11/17
			double theta = acos(2* y / wing.wing_span);
			double numerator = -(2 * N + 1)*cos((2 * N + 1)*theta);
			double denominator = sqrt(pow(wing.wing_span / 2, 2) - pow(y, 2));
			return numerator / denominator;
		};

		// Lets make a quadrature since we have singular endpoints on our integrand, precluding adaptives.
		const int quad_points = 80;
		std::array<double, quad_points> points, weights;
		Quad::gauss_legendre(points, weights);
		for (int idx = 0; idx < quad_points; idx++) {
			Quad::linear_remap(points[idx], weights[idx], -1., 1., -wing.wing_span / 2, wing.wing_span / 2);
		}

		// Compute integ_diff_matrix;
		for (int i = 0; i < number_of_terms; i++) {
			auto ext_coeff = d_3(wing.wing_span*cos(m_collocation_points[i])) 
				/ (2 * Constants::pi() * omega * Constants::i());
			for (int j = 0; j < number_of_terms; j++) {
				// Compute integal.
				integ_diff_matrix(i, j) = ext_coeff * Quad::static_integrate(
					[&](double eta)->std::complex<double> {return sin_deriv_functor(eta, j) * K(wing.wing_span*cos(m_collocation_points[i]) - eta); },
					points, weights, quad_points);
			}
		}

		// Get the LHS matrix:
		LHS_matrix = gamma_matrix - integ_diff_matrix;
		
		// Generate RHS vector:
		if (j == 3) {
			for (int i = 0; i < number_of_terms; i++) {
				double y = wing.wing_span * cos(m_collocation_points[i]);
				RHS_vector[i] = d_3(y);
			}
		}
		else if (j == 5) {
			for (int i = 0; i < number_of_terms; i++) {
				double y = wing.wing_span * cos(m_collocation_points[i]);
				RHS_vector[i] = d_5(y) - U * d_3(y) / (Constants::i() * omega);
			}
		}

		// YAY. We now have everything we need to compute our solution.
		m_solution = LHS_matrix.lu().solve(RHS_vector);
		return;
	}

	std::complex<double> Sclavounos1987::get_solution_vorticity(double y)
	{
		assert(m_solution.size() == number_of_terms);
		std::complex<double> sum(0, 0);
		double theta;
		for (int idx = 0; idx < number_of_terms; idx++) {
			theta = acos(y / wing.wing_span);
			sum += sin((2 * idx + 1) * theta) * m_solution(idx);
		}

		return sum;
	}

	std::complex<double> Sclavounos1987::P(double y)
	{
		double term_1, term_2;
		const int n_points = 10;

		auto integrand_1 = [y](const double t) -> double{
			return exp(-y*t) * (sqrt(t*t - 1) - t) / t;
		};
		auto integrand_2 = [y](const double t) -> double{
			return exp(-y*t) * (sqrt(1 - t*t) - 1) / t;
		};
		
		const int n_pts = 10;
		std::array<double, n_pts> points1, weights1, points2, weights2;
		Quad::gauss_legendre(points1, weights1);
		points2 = points1;
		weights2 = weights1;
		for (int idx = 0; idx < n_pts; idx++) {
			Quad::exponential_remap(points1[idx], weights1[idx], 1.);
			Quad::linear_remap(points2[idx], weights2[idx], -1., 1., 0., 1.);
		}
		// INT term 1 between 1, infty
		term_1 = Quad::static_integrate(integrand_1, points1, weights1, n_points);
		// INT term 2 between 0, 1
		term_2 = Quad::static_integrate(integrand_2, points2, weights2, n_points);

		return term_1 + Constants::i() * term_2;
	}

	double Sclavounos1987::Ei(double y)
	{
		double integral;
		auto integrand = [](double t) -> double {
			return exp(-t) / t;
		};
		
		const int n_pts = 10;
		std::array<double, n_pts> points, weights;
		Quad::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			Quad::exponential_remap(points[idx], weights[idx], -y);
		}
		integral = Quad::static_integrate(integrand, points, weights, n_pts);

		return -1 * integral;
	}

}
