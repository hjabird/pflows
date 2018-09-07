#include "VortexParticleKernels.h"

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
