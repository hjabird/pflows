#pragma once
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.h

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

#include <functional>
#include <complex>
#include <vector>

#include <Eigen/Dense>

#include "WingProjectionGeometry.h"

namespace mFlow {
	class Sclavounos1987
	{
	public:
		Sclavounos1987();
		~Sclavounos1987();

		double U; // Free stream velocity.
		double omega; // Perturbation frequency.
		int j; // Analysis type. 3 or 5.
		WingProjectionGeometry wing; // Geometry definition

		int number_of_terms; // Number of terms in fourier sine series.
		void compute_solution(); // Compute solution.
		std::complex<double> get_solution_vorticity(double y); // Once a solution has been computed values can be retrieved.
		std::complex<double> get_solution_vorticity_deriv(double y); // Rate of change of vorticity wrt/ y.
		std::complex<double> compute_lift_coeff(double heave_added_mass); // Solution must be computed first.


		std::complex<double> F(double y); // Eq5.2 - needs solution first.

		std::vector<double> m_collocation_points; // In 0 to pi (in terms of theta

		double get_elliptic_added_mass_coefficient(double a, double b); // Get added mass for ellipse

	protected:

		Eigen::Matrix<std::complex<double>, -1, 1> m_solution;
		Eigen::Matrix<std::complex<double>, -1, -1> m_gammaprime_K_matrix;

		// Setup
		void compute_collocation_points();

		// See paper for function definitions. Maths.
		std::complex<double> d_3(double y); // Eq4.3
		std::complex<double> d_5(double y); // Eq4.8
		std::complex<double> P(double y); // Eq3.21

		std::complex<double> K(double y); // Eq3.20
		double K_term1(double y); // Singular part of K including numerator term
		double K_term1_singularity(double y); // Just the 1 / (sing_pos - x) term
		double K_term1_singularity_integral(double singularity_pos); // The integral of the singular bit over the span.
		double K_term1_numerator(double y); // The part of the singular term that is nice.
		std::complex<double> K_term2(double y); // The E1 part of K
		std::complex<double> K_term3(double y); // The E_1 and P parts of K

		std::complex<double> integrate_gammaprime_K(double y, int k);
		std::complex<double> integrate_gammaprime_K_term1(double y, int k); // 1/(y-eta) term
		std::complex<double> integrate_gammaprime_K_term2(double y, int k); // E_1 term
		std::complex<double> integrate_gammaprime_K_term3(double y, int k); // P term

		double dtheta_dy(double y); // derivative of theta with respect to y.
		double dsintheta_dy(double y, int k); // derivative of sin((2k + 1)theta) wrt/ y.
		double dsintheta_dtheta(double theta, int k); // derivative of sin((2k+1)theta) wrt/ theta

		// Get a quadrature comprised of two polynomial regions split at split theta.
		// Quadrature is for 0 -> pi. Outputs points_lower, weights_lower, points_upper, weights_upper
		std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
			get_split_quad(int total_pts, double split_theta);
	};
}
