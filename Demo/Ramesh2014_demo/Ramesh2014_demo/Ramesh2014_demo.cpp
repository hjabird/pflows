// Ramesh2014.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <string>

#include <HBTK/AerofoilGenerators.h>
#include <HBTK/AerofoilGeometry.h>
#include <HBTK/AerofoilParser.h>
#include <HBTK/CartesianPlane.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/CsvWriter.h>
#include <HBTK/DoubleTable.h>
#include <HBTK/Generators.h>
#include <HBTK/Paths.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkUnstructuredMeshHolder.h>
#include <HBTK/VtkWriter.h>

#include "Ramesh2014.h"
#include "CanonicalFuntions.h"

int main()
{
	std::cout << "Ramesh2014 Demo. (c) HJAB\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";
	std::cout << "Current working directory: " << HBTK::Paths::current_working_directory() << "\n";
	mFlow::Ramesh2014 sim;

	HBTK::AerofoilGeometry foil;// = HBTK::AerofoilGenerators::sd7003();
	HBTK::CubicSpline1D camber = foil.get_camber_spline();

	sim.semichord = 0.05;
	sim.pitch_location = -1;
	sim.delta_t = 0.00375;
	sim.free_stream_velocity = HBTK::CartesianVector2D({ 0.4, 0 });
	sim.external_purturbation = [](HBTK::CartesianPoint2D p, double t) { 
		//return HBTK::CartesianVector2D({0, sin((p.x()*4 - 2*t)*HBTK::Constants::pi())}); 
		return HBTK::CartesianVector2D({ 0, 0 });
	};
	std::unique_ptr<mFlow::CanonicalFunction> ramp_AoA, ramp_Z, mean_AoA, sum;
	ramp_AoA = std::make_unique<mFlow::Harmonic>(6 * 0.393, -0*0.09477, 0);
	//ramp_AoA = std::make_unique<mFlow::EldredgeSmoothRamp>(1, 3, 4, 6, 11, 0.7854, 1, sim.free_stream_velocity.magnitude());
	mean_AoA = std::make_unique<mFlow::ConstantValue>(HBTK::Constants::degrees_to_radians(4.));
	ramp_Z = std::make_unique<mFlow::Harmonic>(3.2, 0.025, 0);
	sum = std::make_unique<mFlow::CanonicalFunctionSum>(std::move(ramp_AoA), std::move(mean_AoA));
	sim.foil_AoA = [&](double t) { return sum->f(t); };// 0.5 * sin(0.1 * 2 * 3.141 * t); };
	sim.foil_dAoAdt = [&](double t) { return sum->dfdx(t); };// 0.5 * 0.1 * 2 * 3.141*cos(0.1 * 2 * 3.141*t); };
	sim.foil_Z = [&](double t) { return ramp_Z->f(t); };// 0.05 * sin(0.1 * 2 * 3.141 * t); };
	sim.foil_dZdt = [&](double t) { return ramp_Z->dfdx(t); };//0.05 * 0.1*2*3.141*cos(0.1*2*3.141*t); };
	sim.camber_line = [&](double x) { return camber((x + 1) / 2); };
	sim.camber_slope = [&](double x) { return camber.derivative((x + 1) / 2);  };

	sim.number_of_fourier_terms = 8;
	sim.wake_self_convection = true;
	sim.lev_shedding = true;
	sim.critical_leading_edge_suction = 0.08;
	try {
		sim.initialise();
	}
	catch (std::domain_error & e) {
		std::cout << "ERROR:\n";
		std::cout << e.what() << "\n";
		return 1;
	}

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

	const int n_steps = 600;
	for (int i = 0; i < n_steps; i++) {
		std::cout << "\rStep " << i + 1 << " of " << n_steps << "        ";
		sim.advance_one_step();

		step_data["Step"].emplace_back(i + 1);
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

		// Particles
		if (i % 100 == 0) {
			std::string te_particle_name = "output/te_particle_step_" + std::to_string(i) + ".vtu";
			std::ofstream te_particles_file(te_particle_name, std::ios::binary);
			std::string le_particle_name = "output/le_particle_step_" + std::to_string(i) + ".vtu";
			std::ofstream le_particles_file(le_particle_name, std::ios::binary);
			std::string foil_name = "output/foil_step_" + std::to_string(i) + ".vtu";
			std::ofstream foil_file(foil_name, std::ios::binary);
			HBTK::Vtk::VtkUnstructuredDataset data_tev, data_foil;
			int npts = sim.number_of_te_particles();
			HBTK::CartesianPlane plane(
				HBTK::CartesianPoint3D({ 0,0,0 }),
				HBTK::CartesianPoint3D({ 1,0,0 }),
				HBTK::CartesianPoint3D({ 0,0,1 })
			);
			sim.m_te_vortex_particles.save_to_vtk(te_particles_file, plane);
			sim.m_le_vortex_particles.save_to_vtk(le_particles_file, plane);
			auto positions = HBTK::linspace(-1, 1, 30);
			for (int i = 0; i < 30; i++) {
				HBTK::CartesianPoint2D pnt = sim.foil_coordinate(positions[i]);
				data_foil.mesh.points.push_back(plane(pnt));
			}
			for (int i = 0; i < (int)data_foil.mesh.points.size() - 1; i++) {
				data_foil.mesh.cells.push_back({ 3, std::vector<int>({i, i + 1}) });
			}
			HBTK::Vtk::VtkWriter writer;
			writer.appended = false;
			writer.open_file(foil_file, HBTK::Vtk::VtkWriter::UnstructuredGrid);
			writer.write_piece(foil_file, data_foil);
			writer.close_file(foil_file);
		}
	}
	HBTK::CsvWriter csv_writer;
	csv_writer.precision = 6;
	csv_writer.string_limiter = "";
	try {
		csv_writer.write("output/step_data.csv", step_data);
	}
	catch (std::exception& e) {
		std::cout << e.what() << "\n";
	}
	return 0;
}

