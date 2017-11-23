#include "stdafx.h"
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
#include "../../HBTK/HBTK/Constants.h"
#include "../../HBTK/HBTK/Checks.h"

namespace mFlow {
	namespace Common {

		std::complex<double> Common::Hankle2_0(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			return boost::math::cyl_hankel_2(0, k_l);
		}

		std::complex<double> Common::Hankle2_1(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			return boost::math::cyl_hankel_2(1, k_l);
		}

		std::complex<double> Common::Bessel_0(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			return boost::math::cyl_bessel_j(k_l, 0);
		}

		std::complex<double> Common::Bessel_1(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			return boost::math::cyl_bessel_j(k_l, 1);
		}

		double Exponential_int_Ei(double t) {
			assert(HBTK::check_valid(t));
			return boost::math::expint(t);
		}


		std::complex<double> Common::Theodorsen_function(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			return (k_l != 0 ? (Hankle2_1(k_l)) / (Hankle2_1(k_l) + HBTK::Constants::i() * Hankle2_0(k_l)) : 1.);
		}

		std::complex<double> Common::Sears_function(double k_l)
		{
			assert(HBTK::check_valid(k_l));
			std::complex<double> term_11, term_12, term_2;
			term_11 = Theodorsen_function(k_l);
			term_12 = Bessel_0(k_l) + HBTK::Constants::i() * Bessel_1(k_l);
			term_2 = HBTK::Constants::i() * Bessel_1(k_l);
			return term_11 * term_12 + term_2;
		}
	}
}
