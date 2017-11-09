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
		WingProjectionGeometry wing; // Geometry definition (really only chord is used.)

		int number_of_terms; // Number of terms in fourier sine series.
		void compute_solution(); // Compute solution.
		std::complex<double> get_solution_vorticity(double y); // Once a solution has been computed values can be retrieved.

		// Helper
		std::complex<double> P(double y); // Eq3.21
		double Ei(double y); // Used in Eq3.20, https://en.wikipedia.org/wiki/Exponential_integral


	protected:

		Eigen::VectorXd m_fourier_coefficients;
		Eigen::VectorXd m_collocation_points;
		Eigen::Matrix<std::complex<double>, -1, 1> m_solution;

		// Setup
		void compute_collocation_points();

		// See paper for function definitions. Maths.
		std::complex<double> d_3(double y); // Eq4.3
		std::complex<double> d_5(double y); // Eq4.8
		std::complex<double> K(double y); // Eq3.20

	};
}
