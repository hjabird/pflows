#pragma once
/*////////////////////////////////////////////////////////////////////////////
PrandtlLiftingLine.h

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

#include <functional>

#include "WingProjectionGeometry.h"

namespace mFlow {
	class PrandtlLiftingLine {
	public:
		WingProjectionGeometry wing;

		// Incidence with respect to Y_Global
		std::function<double(double)> incidence;
		std::function<double(double)> zero_lift_incidence;

		int n_terms;

		double lift_coefficient();

	private:
		double y_to_theta(double y) const;
		double theta_to_y(double theta);

		std::vector<double> generate_collocation_points();

		double integral(int k, double theta) const;
	};
}