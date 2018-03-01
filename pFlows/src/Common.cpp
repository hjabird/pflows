#include "Common.h"
/*////////////////////////////////////////////////////////////////////////////
Common.cpp

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

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <HBTK/Constants.h>
#include <HBTK/Checks.h>

namespace mFlow {
	namespace Common {

		std::complex<double> hankel2_0(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			return boost::math::cyl_hankel_2(0, k_l);
		}

		std::complex<double> hankel2_1(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			return boost::math::cyl_hankel_2(1, k_l);
		}

		std::complex<double> bessel_0(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			return boost::math::cyl_bessel_j(k_l, 0);
		}

		std::complex<double> bessel_1(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			return boost::math::cyl_bessel_j(k_l, 1);
		}

		double exponential_int_Ei(double t) {
			assert(HBTK::check_finite(t));
			assert(t != 0);
			return boost::math::expint(t);
		}

		double exponential_int_E1(double t) {
			assert(HBTK::check_finite(t));
			return boost::math::expint(1, t);
		}

		std::complex<double> theodorsen_function(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			return (k_l != 0 ? (hankel2_1(k_l)) / (hankel2_1(k_l) + HBTK::Constants::i() * hankel2_0(k_l)) : 1.);
		}

		std::complex<double> sears_function(double k_l)
		{
			assert(HBTK::check_finite(k_l));
			std::complex<double> term_11, term_12, term_2;
			term_11 = theodorsen_function(k_l);
			term_12 = bessel_0(k_l) + HBTK::Constants::i() * bessel_1(k_l);
			term_2 = HBTK::Constants::i() * bessel_1(k_l);
			return term_11 * term_12 + term_2;
		}
	}
}
