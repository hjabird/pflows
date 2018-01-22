#include "stdafx.h"
#include "Guermond1991.h"
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.cpp

Equations from the paper "A unified unsteady lifting-line theory",
J.-L. Guermond and A. Sellier, 1991, J. Fluid. Mech.

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

#include <HBTK/GaussLegendre.h>
#include <HBTK/Integrators.h>
#include <HBTK/Remaps.h>
#include <HBTK/Constants.h>

#include "Common.h"

namespace mFlow {

	double Guermond1991::k_l_function(double k, double c)
	{
		return 0.5 * k * c;
	}

	std::complex<double> Guermond1991::L0_function(std::function<std::complex<double>(double)>& f,
		double k_l, double c_l, double c_t)
	{
		std::complex<double> term_11, term_12;
		std::complex<double> term_21, term_22;

		term_11 = -4.0 * Common::theodorsen_function(k_l);

		const int n_p = 10;
		std::array<double, n_p> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int i = 0; i < n_p; i++) {
			HBTK::telles_cubic_remap(points[i], weights[i], 1.0);
			HBTK::linear_remap(points[i], weights[i], -1., 1., c_l, c_t);
		}

		auto L0_i1 = [&](double xi) {
			return L0_integrand_1(f, c_l, c_t, xi);
		};
		term_12 = HBTK::static_integrate(L0_i1, points, weights, n_p);

		term_21 = (8.0 * HBTK::Constants::i() * k_l) / (c_t - c_l);

		auto L0_i2 = [&](double xi) {
			return L0_integrand_2(f, c_l, c_t, xi);
		};
		term_22 = HBTK::static_integrate(L0_i2, points, weights, n_p);

		return term_11 * term_12 - term_21 * term_22;
	}

	std::complex<double> Guermond1991::M0_function(std::function<std::complex<double>(double)>& f,
		double k_l, double c_l, double c_t)
	{
		std::complex<double> term_11, term_12, term_21, term_22, term_31, term_32;

		term_11 = -(c_t - c_l) * Common::theodorsen_function(k_l);

		auto M0_i1 = [&](double xi) {
			return M0_integrand_1(f, c_l, c_t, xi);
		};

		const int n_p = 10;
		std::array<double, n_p> points, weights;
		HBTK::gauss_legendre(points, weights);
		std::array<double, n_p> points_sato(points), weights_sato(weights);
		for (int i = 0; i < n_p; i++) {
			HBTK::linear_remap(points[i], weights[i], -1., 1., c_l, c_t);
			HBTK::sato_remap<3>(points_sato[i], weights_sato[i], 1.);
			HBTK::linear_remap(points_sato[i], weights_sato[i], -1., 1., c_l, c_t);
		}
		term_12 = HBTK::static_integrate(M0_i1, points_sato, weights_sato, n_p);

		term_21 = c_t - c_l;

		auto M0_i2 = [&](double xi) {
			return M0_integrand_2(f, c_l, c_t, xi);
		};
		term_22 = HBTK::static_integrate(M0_i2, points, weights, n_p);

		term_31 = HBTK::Constants::i() * k_l * (c_t - c_l);

		auto M0_i3 = [&](double xi) {
			return M0_integrand_3(f, c_l, c_t, xi);
		};
		term_32 = HBTK::static_integrate(M0_i3, points, weights, n_p);

		return term_11 * term_12 + term_21 * term_22 + term_31 * term_32;
	}

	std::complex<double> Guermond1991::l_1s_function(double w_1s, double c, double k_l, double K)
	{
		const std::complex<double> i(0, 1);
		std::complex<double> term_11, term_12, term_13, term_14;

		term_11 = 2 * HBTK::Constants::pi() * c;
		term_12 = exp(-i * k_l * (1 - 2 * K));
		term_13 = Common::sears_function(k_l);
		term_14 = w_1s;
		return term_11 * term_12 * term_13 * term_14;
	}

	std::complex<double> Guermond1991::m_1s_function(double w_1s, double c, double k_l, double K)
	{
		std::complex<double> term_11, term_12, term_13, term_14;

		term_11 = HBTK::Constants::pi() *(c*c) / 2.0;
		term_12 = exp(-HBTK::Constants::i() * k_l * (1 - 2 * K));
		term_13 = Common::sears_function(k_l);
		term_14 = w_1s;
		return term_11 * term_12 * term_13 * term_14;
	}

	std::complex<double> Guermond1991::l_1r_function(double c, double k_l, double K, double r_curv,
		double sin_lifting_line_angle, double G, double dGdy, double l0, double dl0dy,
		double l0_prime, double dl0_primedy)
	{
		std::complex<double> term_1, term_21, term_22, term_23,
			term_211, term_212, term_213, term_221, term_222, term_231, term_232;

		term_1 = 2.0 * HBTK::Constants::pi() * c;

		term_211 = exp(-HBTK::Constants::i() * k_l * (1 - 2 * K)) * Common::sears_function(k_l);
		term_212 = log(c) - log(k_l) - log(2.0) - HBTK::Constants::euler() - 0.5 * HBTK::Constants::i() * HBTK::Constants::pi();
		term_213 = curly_A_function(G, dGdy, sin_lifting_line_angle, r_curv);
		term_21 = term_211 * term_212 * term_213;

		term_221 = (Common::theodorsen_function(k_l) - 1.0) / (HBTK::Constants::i() * k_l) + 1. - 2 * K;
		term_222 = curly_A_function(l0, dl0dy, sin_lifting_line_angle, r_curv);
		term_22 = term_221 * term_222;

		term_231 = 2 / c;
		term_232 = curly_A_function(l0_prime, dl0_primedy, sin_lifting_line_angle, r_curv);
		term_23 = term_231 * term_232;

		return term_1 * (term_21 + term_22 - term_23);
	}

	std::complex<double> Guermond1991::m_1r_function(double c, double k_l, double K, double r_curv,
		double sin_lifting_line_angle, double G, double dGdy, double l0, double dl0dy, double l0_prime,
		double dl0_primedy, double l0_doubleprime, double dl0_doubleprimedy)
	{
		std::complex<double> term_1,
			term_21, term_22, term_23, term_24,
			term_211, term_212, term_213,
			term_221, term_222,
			term_231, term_232,
			term_241, term_242;

		term_1 = 0.5 * HBTK::Constants::pi() * c * c;

		term_211 = exp(-HBTK::Constants::i() * k_l * (1 - 2 * K)) * Common::sears_function(k_l);
		term_212 = log(c) - log(k_l) - log(2.0) - HBTK::Constants::euler() - 0.5 * HBTK::Constants::i() * HBTK::Constants::pi();
		term_213 = curly_A_function(G, dGdy, sin_lifting_line_angle, r_curv);
		term_21 = term_211 * term_212 * term_213;

		term_221 = 0.5 * (Common::theodorsen_function(k_l) - 1.0) / (HBTK::Constants::i() * k_l) + pow(1. - 2 * K, 2) - 0.5;
		term_222 = curly_A_function(l0, dl0dy, sin_lifting_line_angle, r_curv);
		term_22 = term_221 * term_222;

		term_231 = 2 * (2 * K - 1.) / c;
		term_232 = curly_A_function(l0_prime, dl0_primedy, sin_lifting_line_angle, r_curv);
		term_23 = term_231 * term_232;

		term_241 = 2. / (c*c);
		term_242 = curly_A_function(l0_doubleprime, dl0_doubleprimedy, sin_lifting_line_angle, r_curv);
		term_24 = term_241 * term_242;

		return term_1 * (term_21 + term_22 + term_23 + term_24);
	}

	std::complex<double> Guermond1991::total_lift(double cos_lifting_line_angle, std::complex<double> l0, std::complex<double> l_1s, std::complex<double> l_1r, double A)
	{
		return cos_lifting_line_angle * (l0 + (l_1s + l_1r) / A);
	}

	std::complex<double> Guermond1991::total_moment(double cos_lifting_line_angle, std::complex<double> m0, std::complex<double> m_1s, std::complex<double> m_1r, double A)
	{
		return cos_lifting_line_angle * (m0 + (m_1s + m_1r) / A);
	}

	std::complex<double> Guermond1991::curly_A_function(std::complex<double> f, std::complex<double> dfdy,
		double sin_lifting_line_angle, double r_curv)
	{
		std::complex<double> term_1, term_2;

		term_1 = f / (4 * HBTK::Constants::pi() * r_curv);
		term_2 = dfdy * sin_lifting_line_angle / (2 * HBTK::Constants::pi());

		return term_1 + term_2;
	}

	std::complex<double> Guermond1991::L0_integrand_1(std::function<std::complex<double>(double)> f, double c_l,
		double c_t, double xi)
	{
		return sqrt((xi - c_l) / (c_t - xi)) * f(xi);
	}

	std::complex<double> Guermond1991::L0_integrand_2(std::function<std::complex<double>(double)> f, double c_l,
		double c_t, double xi)
	{
		return sqrt((xi - c_l) * (c_t - xi)) * f(xi);
	}

	std::complex<double> Guermond1991::M0_integrand_1(std::function<std::complex<double>(double)> f, double c_l,
		double c_t, double xi)
	{
		return L0_integrand_1(f, c_l, c_t, xi);
	}

	std::complex<double> Guermond1991::M0_integrand_2(std::function<std::complex<double>(double)> f, double c_l,
		double c_t, double xi)
	{
		double term_11, term_121, term_122;

		term_11 = sqrt((xi - c_l) / (c_t - xi));
		term_121 = 4.0 / (c_t - c_l);
		term_122 = sqrt((xi - c_l) * (c_t - xi));

		return (term_11 - term_121 * term_122) * f(xi);
	}

	std::complex<double> Guermond1991::M0_integrand_3(std::function<std::complex<double>(double)> f, double c_l,
		double c_t, double xi)
	{
		std::complex<double> term_1;
		double term_21;
		double term_221;
		double term_222;

		const int n_p = 10;
		std::array<double, n_p> points, weights;
		HBTK::gauss_legendre(points, weights);
		/*Quad::apply_remap(points, weights, Quad::linear_remap, -1.0, 1.0, c_l, xi);*/
		for (int i = 0; i < n_p; i++) {
			HBTK::linear_remap(points[i], weights[i], -1., 1., c_l, xi);
		}
		term_1 = HBTK::static_integrate(f, points, weights, n_p);

		term_21 = 1.0 / sqrt((xi - c_l) * (c_t - xi));
		term_221 = 8.0 / ((c_t - c_l)*(c_t - c_l));
		term_222 = sqrt((xi - c_l) * (c_t - xi));

		return term_1 * (term_21 - term_221*term_222);
	}

}