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
#include "../../HBTK/HBTK/Checks.h"

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
		// After Eq 7.1 states these as local maxima of sin((N+1)theta). This doesn't work...
		m_collocation_points.resize(number_of_terms);
		const double hpi = HBTK::Constants::pi() / 2;
		for (int idx = 0; idx < number_of_terms; idx++) {
			m_collocation_points[idx] = (hpi + idx * HBTK::Constants::pi()) / (2 * number_of_terms);
			assert(m_collocation_points[idx] <= HBTK::Constants::pi() && m_collocation_points[idx] >= 0.);
		}
		return;
	}

	std::complex<double> Sclavounos1987::d_3(double y)
	{
		assert(y <= abs(wing.semispan()));
		double l; // semichord
		double v; // Normalised frequency

		l = wing.semichord(y);
		v = omega / U;

		if (omega == 0) { return -2 * HBTK::Constants::pi() * l; }

		auto numerator = 4 * U * exp(-HBTK::Constants::i() * l * v).real(); 
		auto denominator = HBTK::Constants::i() * Common::Hankle2_0(v * l) + Common::Hankle2_1(v * l);
		assert(HBTK::check_finite(numerator));
		assert(HBTK::check_finite(denominator));
		return numerator / denominator;
	}

	std::complex<double> Sclavounos1987::d_5(double y)
	{
		return -0.5 * wing.semichord(y) * d_3(y);
	}

	std::complex<double> Sclavounos1987::K(double y)
	{
		assert(y != 0.);

		if (omega == 0.) { return 1. / (2. * y); } // Eq 5.4

		std::complex<double> term_11, term_12, 
			term_121, term_122, term_123;
		auto v = omega / U;

		term_11 = 0.5 * (y >= 0 ? 1. : -1.);

		term_121 = exp(-v * abs(y)) / abs(y);
		term_122 = HBTK::Constants::i() * v * Common::Exponential_int_Ei(v * abs(y));
		term_123 = v * P(v * abs(y));
		term_12 = term_121 - term_122 + term_123;
		
		assert(HBTK::check_finite(term_11));
		assert(HBTK::check_finite(term_12));
		return term_11 * term_12;
	}

	std::complex<double> Sclavounos1987::F(double y)
	{
		assert(y <= abs(wing.semispan()));

		const int n_pts = 20;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (auto i = 0; i < n_pts; i++) {
			HBTK::linear_remap(points[i], weights[i], -1., 1., -wing.semispan(), wing.semispan());
		}

		auto F_res = -1. * HBTK::static_integrate( 
			[&](double eta){
			return get_solution_vorticity_deriv(eta) * K(y - eta); },	
			points, weights, n_pts) / (wing.span * 2 * HBTK::Constants::pi() * HBTK::Constants::i() * omega);

		assert(HBTK::check_finite(F_res));
		return F_res;
	}


	double Sclavounos1987::dtheta_dy(double y, int N)
	{
		assert(N >= 0);
		assert( y <= abs(wing.semispan()) );
		auto result = -1. / sqrt(pow(wing.semispan(), 2) - pow(y, 2));
		assert(HBTK::check_finite(result));
		return result;
	}


	double Sclavounos1987::dfsintheta_dy(double y, int N)
	{
		assert(y <= abs(wing.semispan()));
		
		auto theta = acos(y / wing.semispan());
		auto dtdy = dtheta_dy(y, N);
		auto dGammadt = (2 * N + 1) * cos((2 * N + 1) * theta);

		return dGammadt * dtdy;
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
		// For the probem in for Ax = b, solve for x, this is A.
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
				gamma_matrix(i, j) = sin((2 * j + 1) * m_collocation_points[i]);
			}
		}

		// Lets make a quadrature since we have singular endpoints on our integrand, precluding adaptives.
		const int quad_points = 80;
		std::array<double, quad_points> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < quad_points; idx++) {
			//HBTK::telles_cubic_remap(points[idx], weights[idx], 0.);
			HBTK::linear_remap(points[idx], weights[idx], -1., 1., -wing.semispan(), wing.semispan());
		}

		// Compute integ_diff_matrix;
		auto one_over_integration_interval = 1 / wing.span;
		for (int i = 0; i < number_of_terms; i++) {

			auto theta = m_collocation_points[i];
			auto y_position = wing.semispan()*cos(theta);
			std::complex<double> ext_coeff;
			if (omega != 0) { 
				ext_coeff = one_over_integration_interval * d_3(y_position)
					/ (2 * HBTK::Constants::pi() * omega * HBTK::Constants::i());
			}
			else {
				ext_coeff = wing.semichord(y_position);
			}

			for (int j = 0; j < number_of_terms; j++) 
			{
				auto integral = HBTK::static_integrate(
					[&](double eta)->std::complex<double> {
						return dfsintheta_dy(y_position, j) * K(y_position - eta); 
					} , // End lambda
					points, weights, quad_points);

				integ_diff_matrix(i, j) = ext_coeff * integral;
			} // End For in j
		} // End for in i

		LHS_matrix = gamma_matrix - integ_diff_matrix;
		
		// Generate RHS vector:
		if (j == 3) {
			for (int i = 0; i < number_of_terms; i++) {
				auto y = wing.semispan()*cos(m_collocation_points[i]);
				RHS_vector[i] = d_3(y);
			}
		}
		else if (j == 5) {
			for (int i = 0; i < number_of_terms; i++) {
				auto y = wing.semispan()*cos(m_collocation_points[i]);
				RHS_vector[i] = d_5(y) - U * d_3(y) / (HBTK::Constants::i() * omega);
			}
		}

		// YAY. We now have everything we need to compute our solution.
		assert(abs(LHS_matrix.determinant()) >= 1e-5);
		m_solution = LHS_matrix.lu().solve(RHS_vector);
		return;
	}

	std::complex<double> Sclavounos1987::get_solution_vorticity(double y)
	{
		assert(m_solution.size() == number_of_terms);
		assert(y <= abs(wing.semispan()) );

		std::complex<double> sum(0, 0);
		double theta;
		for (int idx = 0; idx < number_of_terms; idx++) {
			theta = acos(y / wing.semispan());
			assert(HBTK::check_finite(m_solution(idx)));
			sum += sin((2 * idx + 1) * theta) * m_solution(idx);
		}

		assert(HBTK::check_finite(sum));
		return sum;
	}

	std::complex<double> Sclavounos1987::get_solution_vorticity_deriv(double y)
	{
		assert(m_solution.size() == number_of_terms);
		assert(y <= abs(wing.semispan()));

		std::complex<double> sum(0, 0);
		double theta;
		for (int idx = 0; idx < number_of_terms; idx++) {
			theta = acos( y / wing.semispan());
			assert(HBTK::check_finite(m_solution(idx)));
			sum += dfsintheta_dy(y, idx) * m_solution(idx);
		}

		assert(HBTK::check_finite(sum));
		return sum;
	}

	std::complex<double> Sclavounos1987::compute_lift_coeff(double heave_added_mass)
	{
		std::complex<double> term_1, term_2,
			term_11, term_12, term_21, term_22;

		auto integrand = [&](double y) -> std::complex<double> {
			auto semichord = wing.semichord(y);
			auto F_3 = F(y);
			auto C = Common::Theodorsen_function((omega / U) * semichord);
			return C * semichord * F_3;
		};

		auto wing_area = wing.area();

		term_11 = -4. / wing_area;
		const int n_pts = 80;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			HBTK::linear_remap(points[idx], weights[idx], -1., 1., -wing.semispan(), wing.semispan());
		}
		term_12 = HBTK::static_integrate(integrand, points, weights, n_pts);
		term_1 = term_11 * term_12;

		term_21 = HBTK::Constants::i() * (omega / U);
		term_22 = heave_added_mass / wing_area;
		term_2 = term_21 * term_22;

		assert(HBTK::check_finite(term_1));
		assert(HBTK::check_finite(term_2));
		return -term_1 - term_2;
	}

	std::complex<double> Sclavounos1987::P(double y)
	{
		double term_1, term_2;

		auto integrand_1 = [y](const double t) -> double{
			return exp(-y*t) * (sqrt(t*t - 1) - t) / t;
		};
		auto integrand_2 = [y](const double t) -> double{
			return exp(-y*t) * (sqrt(1 - t*t) - 1) / t;
		};
		
		const int n_points = 20;
		std::array<double, n_points> points1, weights1, points2, weights2;
		HBTK::gauss_legendre(points1, weights1);
		points2 = points1;
		weights2 = weights1;
		for (int idx = 0; idx < n_points; idx++) {
			HBTK::exponential_remap(points1[idx], weights1[idx], 1.);
			HBTK::linear_remap(points2[idx], weights2[idx], -1., 1., 0., 1.);
		}

		term_1 = HBTK::static_integrate(integrand_1, points1, weights1, n_points);
		term_2 = HBTK::static_integrate(integrand_2, points2, weights2, n_points);
		assert(HBTK::check_finite(term_1));
		assert(HBTK::check_finite(term_2));

		return term_1 + HBTK::Constants::i() * term_2;
	}

	/*
	double Sclavounos1987::E_1(double t)
	{
		assert(t != 0);

		auto integrand = [](double t) -> double {
			return exp(-t) / t;
		};
		
		const int n_pts = 15;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			HBTK::exponential_remap(points[idx], weights[idx], t);
		}
		auto integral = HBTK::static_integrate(integrand, points, weights, n_pts);

		return integral;
	}
	*/
}
