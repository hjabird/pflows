#pragma once

#include <complex>

namespace mFlow {
	namespace Common {

		std::complex<double> Hankle2_0(double); // H^{(2)}_{0}
		std::complex<double> Hankle2_1(double); // H^{(2)}_{1}
		std::complex<double> Bessel_0(double); // J_0
		std::complex<double> Bessel_1(double); // J_1

		std::complex<double> Theodorsen_function(double k_l); // C(k_1) 
		std::complex<double> Sears_function(double k_l); // S(k_1)

	}
}
