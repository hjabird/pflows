#include "McGowan2011.h"
/*////////////////////////////////////////////////////////////////////////////
McGowan2011.cpp

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

#include <HBTK/Constants.h>

#include "Common.h"

mFlow::McGowan2011::McGowan2011()
{
}


mFlow::McGowan2011::~McGowan2011()
{
}

double mFlow::McGowan2011::reduced_frequency()
{
	return frequency * semichord / free_stream_vel;
}

std::complex<double> mFlow::McGowan2011::plunge_amp()
{
	return plunge_amplitude;
}

std::complex<double> mFlow::McGowan2011::dplunge_dt_amp()
{
	return HBTK::Constants::i() * frequency * plunge_amplitude;
}

std::complex<double> mFlow::McGowan2011::ddplunge_ddt_amp()
{
	return - frequency * frequency * plunge_amplitude;
}

std::complex<double> mFlow::McGowan2011::pitch_amp()
{
	return pitch_amplitude * exp(HBTK::Constants::i() * phase_offset);
}

std::complex<double> mFlow::McGowan2011::dpitch_dt_amp()
{
	return HBTK::Constants::i() * frequency * pitch_amplitude * exp(HBTK::Constants::i() * phase_offset);
}

std::complex<double> mFlow::McGowan2011::ddpitch_ddt_amp()
{
	return -frequency * frequency * pitch_amplitude * exp(HBTK::Constants::i() * phase_offset);
}

std::complex<double> mFlow::McGowan2011::qstat_unsteady()
{
	std::complex<double> term_11, term_12, term_121, term_122, term_123;

	term_11 = 2 * HBTK::Constants::pi();
	
	term_121 = -dplunge_dt_amp() / free_stream_vel;
	term_122 = pitch_amp();
	term_123 = semichord * (0.5 - pitch_location) * dpitch_dt_amp() / free_stream_vel;
	term_12 = term_121 + term_122 + term_123;

	return term_11 * term_12;
}

std::complex<double> mFlow::McGowan2011::theodorsen_unsteady()
{
	std::complex<double> term_1, term_2,
		term_11, term_12, term_21, term_22,
		term_121, term_122, term_123,
		term_221, term_222, term_223;

	term_11 = HBTK::Constants::pi() * semichord;
	term_121 = -ddplunge_ddt_amp() / pow(free_stream_vel, 2);
	term_122 = dpitch_dt_amp() / free_stream_vel;
	term_123 = -semichord * pitch_location * ddpitch_ddt_amp() / pow(free_stream_vel, 2);
	term_12 = term_121 + term_122 + term_123;
	term_1 = term_11 * term_12;

	term_21 = 2. * HBTK::Constants::pi() * Common::theodorsen_function(reduced_frequency());
	term_221 = -dplunge_dt_amp() / free_stream_vel;
	term_222 = pitch_amp();
	term_223 = semichord * (0.5 - pitch_location) * dpitch_dt_amp() / free_stream_vel;
	term_22 = term_221 + term_222 + term_223;
	term_2 = term_21 * term_22;

	return term_1 + term_2;
}
