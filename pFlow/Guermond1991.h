#pragma once
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.h

Equations from the paper "A unified unsteady lifting-line theory",
J.-L. Guermond and A. Sellier, 1991, J. Fluid. Mech.

Copyright 2017 HJA Bird

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/////////////////////////////////////////////////////////////////////////////

#include <complex>
#include <functional>

namespace Guermond1991{

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

	// Special functions wrapper for external libraries.
	std::complex<double> Hankle2_0(double); // H^{(2)}_{0}
	std::complex<double> Hankle2_1(double); // H^{(2)}_{1}
	std::complex<double> Bessel_0(double); // J_0
	std::complex<double> Bessel_1(double); // J_1

	std::complex<double> Theodorsen_function(double k_l); // C(k_1) - Eq66
	std::complex<double> Sears_function(double k_l); // S(k_1) - Eq72

	double k_l_function(double k, double c); // Defd between Eq65-Eq66

	// 2D (Section 6.1)
	std::complex<double> L0_function(std::function<std::complex<double>(double)> & f,
		double k_l, double c_l, double c_t); // Eq65
	std::complex<double> M0_function(std::function<std::complex<double>(double)> & f,
		double k_l, double c_l, double c_t); // Eq67

	// 3D corrections (Section 6.2)
	std::complex<double> l_1s_function(double w_1s,	double c, double k_l, double K); // Eq70
	std::complex<double> m_1s_function(double w_1s,	double c, double k_l, double K); // Eq71
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