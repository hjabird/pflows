#pragma once
/*////////////////////////////////////////////////////////////////////////////
Common.h

Equations that are common to multiple papers.

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

namespace mFlow {
	namespace Common {

		std::complex<double> hankel2_0(double); // H^{(2)}_{0}
		std::complex<double> hankel2_1(double); // H^{(2)}_{1}
		std::complex<double> bessel_0(double);	// J_0
		std::complex<double> bessel_1(double);	// J_1

		double exponential_int_Ei(double t);	// Ei(x)
		double exponential_int_E1(double t);	// E_1(x)

		std::complex<double> theodorsen_function(double k_l);	// C(k_1) 
		std::complex<double> sears_function(double k_l);		// S(k_1)

	}
}
