// Ramesh2014.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <string>

#include <HBTK/AerofoilGeometry.h>
#include <HBTK/AerofoilParser.h>
#include <HBTK/CsvWriter.h>
#include <HBTK/DoubleTable.h>
#include <HBTK/Generators.h>

#include "Ramesh2014.h"

int main()
{

	std::cout << "Ramesh2014 Demo. (c) HJAB\n";
	mFlow::Ramesh2014 sim;

	HBTK::AerofoilGeometry foil;
	HBTK::AerofoilParser foil_parser(foil);
	foil_parser.parse("inputs/SD7003.dat");
	HBTK::CubicSpline1D camber = foil.get_camber_spline();

	sim.semichord = 0.5;
	sim.pitch_location = 0;
	sim.delta_t = 0.01;
	sim.free_stream_velocity = 1;
	sim.foil_AoA = [](double t) { return 0.07; };
	sim.foil_dAoAdt = [](double t) { return 0.0; };
	sim.foil_Z = [](double t) { return 0.05 * sin(7.86*t); };
	sim.foil_dZdt = [](double t) { return 0.05 * 7.86 * cos(7.86*t); };
	sim.camber_line = [&](double x) { return camber((x + 1) / 2);  };
	sim.camber_slope = [&](double x) { return camber.derivative((x + 1) / 2);  };

	sim.number_of_fourier_terms = 8;
	sim.initialise();

	HBTK::DoubleTable step_data;
	step_data.add_column("Time");
	step_data.add_column("Step");
	step_data.add_column("Delta_t");
	step_data.add_column("AoA");
	step_data.add_column("h");
	step_data.add_column("Wake_vorticity");
	step_data.add_column("Cl");
	step_data.add_column("Cd");
	step_data.add_column("A0");
	step_data.add_column("A1");
	step_data.add_column("A2");

	for (int i = 0; i < 1000; i++) {
		sim.advance_one_step();

		step_data["Step"].emplace_back(i+1);
		step_data["Time"].emplace_back(sim.time);
		step_data[2].emplace_back(sim.delta_t);
		step_data["AoA"].emplace_back(sim.foil_AoA(sim.time));
		step_data["h"].emplace_back(sim.foil_Z(sim.time));
		step_data[5].emplace_back(sim.total_shed_vorticity());
		double cl, cd;
		std::tie(cl, cd) = sim.aerofoil_lift_and_drag_coefficients();
		step_data["Cl"].emplace_back(cl);
		step_data["Cd"].emplace_back(cd);
		auto fourier_coeff = sim.get_fourier_terms();
		step_data[8].emplace_back(fourier_coeff[0]);
		step_data[9].emplace_back(fourier_coeff[1]);
		step_data[10].emplace_back(fourier_coeff[2]);

		std::string particle_name = "output/particles_step_" + std::to_string(i) + ".vtk";
		std::string foil_name = "output/foil_step_" + std::to_string(i) + ".vtk";
		std::ofstream particles_file(particle_name, std::ios::binary);
		std::ofstream foil_file(foil_name, std::ios::binary);
		int npts = sim.number_of_particles();
		particles_file << "# vtk DataFile Version 2.0\nFile\nASCII\n";
		particles_file << "DATASET UNSTRUCTURED_GRID\n";
		particles_file << "POINTS " << npts << " float\n";
		for (auto particle : sim.m_vortex_particles) {
			particles_file << particle.x << " 0.0 " << particle.y << "\n";
		}

		particles_file << "\nCELLS " << npts << " " << npts * 2 << "\n";
		for (int j = 0; j < npts; j++) {
			particles_file << "1 " << j << "\n";
		}
		particles_file << "\nCELL_TYPES " << npts << "\n";
		for (int j = 0; j < npts; j++) { particles_file << "1\n"; }
		particles_file << "\nPOINT_DATA " << npts << "\n";
		particles_file << "SCALARS vorticity float 1\n";
		particles_file << "LOOKUP_TABLE default\n";
		for (int j = 0; j < npts; j++) {
			particles_file << sim.m_vortex_particles[j].vorticity << "\n";
		}
		particles_file << "VECTORS velocities float\n";
		for (int j = 0; j < npts; j++) {
			particles_file << sim.m_vortex_particles[j].vx << " 0.0 "
				<< sim.m_vortex_particles[j].vy << "\n";
		}



		foil_file << "# vtk DataFile Version 2.0\nFile\nASCII\n";
		foil_file << "DATASET UNSTRUCTURED_GRID\n";
		foil_file << "POINTS 30 float\n";
		double xtmp, ztmp;
		auto positions = HBTK::linspace(-1, 1, 30);
		for (int i = 0; i < 30; i++) {
			std::tie(xtmp, ztmp) = sim.foil_coordinate(positions[i]);
			foil_file << xtmp << " 0.0 " << ztmp << "\n";
		}
		foil_file << "\nCELLS 29 87\n";
		for (int i = 0; i < 29; i++) { foil_file << "2 " << i << " " << i + 1 << "\n"; }
		foil_file << "\nCELL_TYPES 29\n";
		for (int i = 0; i < 29; i++) { foil_file << "3" << "\n"; }

	}
	HBTK::CsvWriter csv_writer;
	csv_writer.precision = 6;
	csv_writer.string_limiter = "";
	csv_writer.write("output/step_data.csv", step_data);
	std::cout << sim.number_of_particles();




	return 0;
}

