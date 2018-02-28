#pragma once

/*////////////////////////////////////////////////////////////////////////////
Plotting.h

Functions for common plotting tasks.

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
#include <HBTK/GnuPlot.h>

namespace mFlow {
	namespace Plotting{

		// Give something like a complex lift coefficient,  plot it over the period given. (Sinusiodal!)
		void plot_complex_over_period(std::complex<double> plot_variable, double oscillation_period,
 HBTK::GnuPlot & plot,
			std::string line_spec);

	}

}

