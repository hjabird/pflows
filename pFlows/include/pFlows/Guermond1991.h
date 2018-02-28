#pragma once
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.h

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

#include <complex>
#include <functional>

namespace mFlow {
	namespace Guermond1991 {

		// This is based on Guermond and Sellier 1991 - A unified unsteady lifting line 
		// theory. Variables and functions are generally made to resemble the notation 
		// used in the paper.

		// x - streamwise dir (+ve x downstream of wing)
		// y - spanwise dir
		// c(y) - chord of wing at y = c_t - c_l
		// c_t(y) - trailing edge x pos.
		// c_l(y) - leading edge x pos.
		// f(xi) - the imposed normal velocity in x direction defd ~Eq65
		// k - chord reduced frequency defd Eq12
		// k_l - local chord reduced frequency (see k1_function)  
		// K - distance of lifting line from leading edge as fraction of chord c(y).
		// l0 - two dimestionsal lift
		// l0_prime - first moment of 2D acceleration potential
		// l0_doubleprime - second momoent of 2D accl pot.
		// r_curv - the radius of curvature of the lifting line at y

		double k_l_function(double k, double c); // Defd between Eq65-Eq66

		// 2D (Section 6.1)
		std::complex<double> L0_function(std::function<std::complex<double>(double)> & f,
			double k_l, double c_l, double c_t); // Eq65
		std::complex<double> M0_function(std::function<std::complex<double>(double)> & f,
			double k_l, double c_l, double c_t); // Eq67

		// 3D corrections (Section 6.2)
		std::complex<double> l_1s_function(double w_1s, double c, double k_l, double K); // Eq70
		std::complex<double> m_1s_function(double w_1s, double c, double k_l, double K); // Eq71
		std::complex<double> l_1r_function(double c, double k_l, double K, double r_curv,
			double sin_lifting_line_angle, double G, double dGdy,
			double l0, double dl0dy, double l0_prime, double dl0_primedy); // Eq73
		std::complex<double> m_1r_function(double c, double k_l, double K, double r_curv,
			double sin_lifting_line_angle, double G, double dGdy,
			double l0, double dl0dy, double l0_prime, double dl0_primedy,
			double l0_doubleprime, double dl0_doubleprimedy); // Eq76
		std::complex<double> total_lift(double cos_lifting_line_angle, std::complex<double> l0,
			std::complex<double> l_1s, std::complex<double> l_1r, double A); // Eq77
		std::complex<double> total_moment(double cos_lifting_line_angle, std::complex<double> m0,
			std::complex<double> m_1s, std::complex<double> m_1r, double A); // Eq78

		// Convenience functions
		std::complex<double> curly_A_function(std::complex<double> f, std::complex<double> dfdy,
			double sin_lifting_line_angle, double r_curv); // Eq74


		// Integrands in sub expressions.
		// Eq65 integrands
		std::complex<double> L0_integrand_1(std::function<std::complex<double>(double)> f,
			double c_l, double c_t, double xi);
		std::complex<double> L0_integrand_2(std::function<std::complex<double>(double)> f,
			double c_l, double c_t, double xi);
		// Eq67 integrands
		std::complex<double> M0_integrand_1(std::function<std::complex<double>(double)> f,
			double c_l, double c_t, double xi);
		std::complex<double> M0_integrand_2(std::function<std::complex<double>(double)> f,
			double c_l, double c_t, double xi);
		std::complex<double> M0_integrand_3(std::function<std::complex<double>(double)> f,
			double c_l, double c_t, double xi);
	}
}
