
#include <iostream>

#include <HBTK/Paths.h>

#include "Guermond1990.h"
#include "WingGenerators.h"

int main(int argc, char* argv[]) {
	std::cout << "Demo of Guermond 1990.\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path();
	std::cout << "Working path: " << HBTK::Paths::current_working_directory();

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, 1, 4);

	mFlow::Guermond1990 sim;
	sim.wing = wing;
	sim.incidence = [](double x) { return 0.1; };
	sim.zero_lift_incidence = [](double x) { return 0; };
	sim.moment_about_midchord = [&](double x) { return sim.moment_about_midchord_no_camber_slope(x); };

	std::cout << "Circulation at middle of wing is " << sim.Gamma_0(0.0) << "\n";

	return 0;
}