#pragma once

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

