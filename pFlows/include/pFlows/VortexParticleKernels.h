#pragma once
/*////////////////////////////////////////////////////////////////////////////
VortexParticleKernels.h

Vortex particle kernels as defined in D. K. Robertson and G. W. Reich, 
“3-D Vortex Particle Aerodynamic Modelling and Trajectory Optimization of 
Perching Manoeuvres,” in 54th AIAA/ASME/ASCE/AHS/ASC Structures, Structural 
Dynamics, and Materials Conference, 2013.

Copyright 2018 HJA Bird

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

#include <cmath>

#include <HBTK/Constants.h>

namespace mFlow {
	enum VortexParticleRegularisation {
		SINGULAR,
		PLANETARY,
		EXPONENTIAL,
		WINCKELMANS,
		TANH
	};

	struct singular_vortex{
		static double g_function(double rho);
		static double f_function(double rho);
	};
	struct planetary_vortex {
		static double g_function(double rho);
		static double f_function(double rho);
	};
	struct exponential_vortex {
		static double g_function(double rho);
		static double f_function(double rho);
	};
	struct winckelmans_vortex {
		static double g_function(double rho);
		static double f_function(double rho);
	};
	struct tanh_vortex {
		static double g_function(double rho);
		static double f_function(double rho);
	};
}

