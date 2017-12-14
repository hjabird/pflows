#pragma once
/*////////////////////////////////////////////////////////////////////////////
Sclavounos1987.h

Equations from the paper "An unsteady lifting-line theory", P.D. Sclavounos,
1987, Journal of Engineering Mathmatics.
https://doi.org/10.1007/BF00127464

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

		// Geometry definition
		WingProjectionGeometry wing;

		// Parameters for simulation
		double U;		// Free stream velocity.
		double omega;	// Perturbation frequency.
		int j;			// Analysis type. 3 or 5, corresponding to heave and pitch respectively.
		std::vector<double> m_collocation_points;	// In 0 to pi (in terms of theta)
		int number_of_terms;	// Number of sin((2k+1)theta) terms used to approximate the vorticity distribution.

		// Assuming that valid input parameters have been set, compute the solution, leaving results in m_solution.
		void compute_solution();

		// Once the solution has been computed, the lift coefficient may be computed. The heave added mass
		// coefficient, given as a_{33} for heave added mass must be given. NB: a_{33} = 2m_{33} / \rho
		std::complex<double> compute_lift_coeff(double heave_added_mass);

		// Once a solution has been computed, the complex vorticity at any point on the lifting line can be 
		// evaluated.
		std::complex<double> solution_vorticity(double y);

		// Once a solution has been computed, the interaction term, F (see eq5.2) can be evaluated
		// at any point on the lifthing line.
		std::complex<double> F(double y);

		// Returns the added mass coefficient for the wing on the assumption it is elliptic.
		double elliptic_added_mass_coefficient();

	protected:

		// The coefficients A_{2k+1} computed using eq5.2.
		Eigen::Matrix<std::complex<double>, -1, 1> m_solution;
		// The matrix integrate(Gamma'(eta) * K(y-eta) deta) was expensive to evaluate and 
		// so is cached.
		Eigen::Matrix<std::complex<double>, -1, -1> m_gammaprime_K_matrix;

		// Setup
		void compute_collocation_points();

		// Unsteady vortex coefficients:
		std::complex<double> d_3(double y); // Eq4.3 - heave
		std::complex<double> d_5(double y); // Eq4.8 - pitch

		// The Kernal described in eq3.20 in terms of K(y) = ...
		std::complex<double> K(double y);		// Eq3.20
		double K_term1(double y);				// Singular part of K including numerator term
		double K_term1_singularity(double y);	// Just the 1 / (sing_pos - x) term
		double K_term1_numerator(double y);		// The part of the singular term that is nice.
		std::complex<double> K_term2(double y); // The E1 part of K
		std::complex<double> K_term3(double y); // The E_1 and P parts of K

		// A function used in K(y) (eq3.20).
		std::complex<double> P(double y);		// Eq3.21

		// Integrate Int(Gamma'(eta) * K(y - eta) dy) from eq5.2, eq5.3.
		// y is the y position, k part of a frequency term - Gamma = Sum(A_{2k+1} sin((2k+1)theta))
		std::complex<double> integrate_gammaprime_K(double y, int k);
		std::complex<double> integrate_gammaprime_K_term1(double y, int k); // 1/(y-eta) term
		std::complex<double> integrate_gammaprime_K_term2(double y, int k); // E_1 term
		std::complex<double> integrate_gammaprime_K_term3(double y, int k); // P term

		// Derivatives of the Gamma term.
		double dtheta_dy(double y);						// derivative of theta with respect to y.
		double dsintheta_dy(double y, int k);			// derivative of sin((2k + 1)theta) wrt/ y.
		double dsintheta_dtheta(double theta, int k);	// derivative of sin((2k+1)theta) wrt/ theta


		// Get a quadrature comprised of two polynomial regions split at split theta.
		// Quadrature is for 0 -> pi. Outputs points_lower, weights_lower, points_upper, weights_upper
		std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
			split_quad(int total_pts, double split_theta);
	};
}
