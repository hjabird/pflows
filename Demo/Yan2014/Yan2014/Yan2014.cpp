// Yan2014.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cmath>
#include <fstream>

#include <HBTK/VtkLegacyWriter.h>

#include "Yan2014.h"


int main()
{
	std::cout << "Yan2014 Demo. (c) HJAB\n";
	mFlow::Yan2014 sim;

	sim.semichord = 1;
	sim.delta_t = 0.01;
	sim.free_stream_velocity = 1;
	sim.foil_AoA = [](double t) { return 0.05; };
	sim.foil_dAoAdt = [](double t) { return 0.0; };
	sim.foil_Z = [](double t) { return  0.0; };
	sim.foil_dZdt = [](double t) { return  0.0; };
	
	for (int i = 0; i < 400; i++) {
		sim.advance_one_step();

		std::string particle_name = "output/particles_step_" + std::to_string(i) + ".vtk";
		std::string foil_name = "output/foil_step_" + std::to_string(i) + ".vtk";
		std::ofstream particles_file(particle_name, std::ios::binary);
		std::ofstream foil_file(foil_name, std::ios::binary);
		int npts = sim.number_of_particles();
		particles_file << "# vtk DataFile Version 2.0\nFile\nASCII\n";
		particles_file << "DATASET UNSTRUCTURED_GRID\n";
		particles_file << "POINTS " << npts  << " float\n";
		for (auto particle : sim.m_vortex_particles) {
			particles_file << particle.X << " 0.0 " << particle.Z << "\n";
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
			particles_file << sim.m_vortex_particles[j].vX << " 0.0 "
				<< sim.m_vortex_particles[j].vZ << "\n";
		}



		foil_file << "# vtk DataFile Version 2.0\nFile\nASCII\n";
		foil_file << "DATASET UNSTRUCTURED_GRID\n";
		foil_file << "POINTS 2 float\n";
		double lex, lez, tex, tez;
		lex = sim.semichord * cos(sim.foil_AoA(sim.time));
		tex = -lex;
		lez = sim.foil_Z(sim.time) - sim.semichord * sin(sim.foil_AoA(sim.time));
		tez = sim.foil_Z(sim.time) + sim.semichord * sin(sim.foil_AoA(sim.time));
		foil_file << tex << " 0.0 " << tez << "\n";
		foil_file << lex << " 0.0 " << lez << "\n";
		foil_file << "\nCELLS 1 3\n2 0 1\n";
		foil_file << "\nCELL_TYPES 1\n3\n";
	}
	std::cout << sim.number_of_particles();




	return 0;
}

