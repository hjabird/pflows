
#include <iostream>

#include <HBTK/Paths.h>

#include "PrandtlLiftingLine.h"
#include "WingGenerators.h"

int main(int argc, char* argv[]) {
	std::cout << "Demo of Prandtl LLT.\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";
	std::cout << "Working path: " << HBTK::Paths::current_working_directory() << "\n";

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, 1, 4);

	mFlow::PrandtlLiftingLine sim;
	sim.wing = wing;
	sim.incidence = [](double x) { return 0.01745 * 4; };
	sim.zero_lift_incidence = [](double x) { return -0.01745 * 1.8; };
	sim.n_terms = 4;
	std::cout << "Lift coeff is " << sim.lift_coefficient() << "\n";

	return 0;
}