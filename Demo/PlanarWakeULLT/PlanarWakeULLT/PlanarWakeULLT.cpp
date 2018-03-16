// Demo for PlanarWakeULLT

#include <iostream>

#include <HBTK/AerofoilGenerators.h>
#include <HBTK/AerofoilGeometry.h>
#include <HBTK/Paths.h>

#include "PlanarWakeULLT.h"
#include "WingGenerators.h"
#include "WingProjectionGeometry.h"

int main(int argc, char* argv[]) {
	std::cout << "ULLT Demo (C) HJA Bird 2018\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";
	std::cout << "Current working directory: " << HBTK::Paths::current_working_directory() << "\n";
	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, 1, 4);
	HBTK::AerofoilGeometry foil = HBTK::AerofoilGenerators::naca_four_digit("0012");

	mFlow::LUALLT_planar_wake sim;

	int num_inner = 5;
	for (int i = 0; i < num_inner; i++) {
		double y_pos = -wing.semispan() + wing.span * (i + 0.5) / (num_inner + 1);
		P
	}




	return 0;
}