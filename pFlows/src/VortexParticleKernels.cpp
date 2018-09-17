#include "VortexParticleKernels.h"
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

double mFlow::singular_vortex::g_function(double rho) { return 1.0; }
double mFlow::singular_vortex::f_function(double rho) { return 0.0; }

double mFlow::planetary_vortex::g_function(double rho) { return (rho < 1 ? pow(rho, 3) : 1); }
double mFlow::planetary_vortex::f_function(double rho) { return rho < 1. ? 3 : 0; }

double mFlow::exponential_vortex::g_function(double rho) { return 1 - f_function(rho) / 3.; }
double mFlow::exponential_vortex::f_function(double rho) { return 3 * exp(-pow(rho, 3)); }

double mFlow::winckelmans_vortex::g_function(double rho) {
	double num = (rho * rho + 2.5) * pow(rho, 3);
	double den = pow(rho * rho + 1, 2.5);
	return num / den;
}
double mFlow::winckelmans_vortex::f_function(double rho) { return 7.5 / pow(rho * rho + 1, 3.5); }

double mFlow::tanh_vortex::g_function(double rho) { return tanh(rho * rho * rho); }
double mFlow::tanh_vortex::f_function(double rho) { return 3 * pow(1 / cosh(pow(rho, 3)), 2); }
