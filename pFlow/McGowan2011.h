#pragma once
/*////////////////////////////////////////////////////////////////////////////
McGowan2011.h

Equations from the paper "Investigations of Lift-Based Pitch-Plunge Equivalence
for Airfoils at Low Reynolds Numbers", AIAA Journal, 2011.
DOI: 10.2514/1.J050924

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
	class McGowan2011
	{
	public:
		McGowan2011();
		~McGowan2011();

		double frequency;			// omega
		double free_stream_vel;		// V_infty
		double plunge_amplitude;
		double pitch_amplitude;
		double phase_offset;		// In radians
		double pitch_location;		// Between -1, 1 as distance from LE.
		double semichord;

		double reduced_frequency();		// k -> (omega c) / (2 V_infty)

		// Resulting coefficient of lift at t=0 for quasi-steady thin aerofoil theory.
		std::complex<double> qstat_unsteady();
		// Resulting coefficient of lift at t=0 for Theodorsen's theory.
		std::complex<double> theodorsen_unsteady();



	protected:
		// Complex amplitude of pitch and plunge functions. Should be multiplied
		// through by exp(i omega t) and real part is taken to get actual value.
		std::complex<double> plunge_amp();
		std::complex<double> dplunge_dt_amp();
		std::complex<double> ddplunge_ddt_amp();
		std::complex<double> pitch_amp();
		std::complex<double> dpitch_dt_amp();
		std::complex<double> ddpitch_ddt_amp();

	};
}

