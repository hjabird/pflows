#include "stdafx.h"
#include "Common.h"

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "../../HBTK/HBTK/Constants.h"

namespace mFlow {
	namespace Common {

		std::complex<double> Common::Hankle2_0(double k_l)
		{
			return boost::math::cyl_hankel_2(0, k_l);
		}

		std::complex<double> Common::Hankle2_1(double k_l)
		{
			return boost::math::cyl_hankel_2(1, k_l);
		}

		std::complex<double> Common::Bessel_0(double k_l)
		{
			return boost::math::cyl_bessel_j(k_l, 0);
		}

		std::complex<double> Common::Bessel_1(double k_l)
		{
			return boost::math::cyl_bessel_j(k_l, 1);
		}


		std::complex<double> Common::Theodorsen_function(double k_l)
		{
			return (Hankle2_1(k_l)) / (Hankle2_1(k_l) + Constants::i() * Hankle2_0(k_l));
		}

		std::complex<double> Common::Sears_function(double k_l)
		{
			std::complex<double> term_11, term_12, term_2;
			term_11 = Theodorsen_function(k_l);
			term_12 = Bessel_0(k_l) + Constants::i() * Bessel_1(k_l);
			term_2 = Constants::i() * Bessel_1(k_l);
			return term_11 * term_12 + term_2;
		}
	}
}
