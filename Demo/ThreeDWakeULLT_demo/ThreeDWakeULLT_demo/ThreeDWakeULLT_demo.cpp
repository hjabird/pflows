/*////////////////////////////////////////////////////////////////////////////
ThreeDWakeULLT_demo.cpp

A demonstration of the ThreeDWakeULLT method.

Copyright 2017 HJA Bird

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

#include <iostream>
#include <omp.h>

#include <HBTK/AerofoilGenerators.h>
#include <HBTK/AerofoilGeometry.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/Constants.h>
#include <HBTK/CsvWriter.h>
#include <HBTK/DoubleTable.h>
#include <HBTK/Generators.h>
#include <HBTK/Paths.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

#include "ThreeDWakeULLT.h"
#include "Ramesh2014.h"
#include "WingGenerators.h"
#include "WingProjectionGeometry.h"
#include "CanonicalFunctions.h"


int main(int argc, char* argv[]) {
	std::cout << "ULLT Demo (C) HJA Bird 2018\n";
#ifdef _OPENMP
	std::cout << "Compiled with OpenMP (Max threads = " << omp_get_max_threads() << ").\n";
#endif
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";
	std::cout << "Current working directory: " << HBTK::Paths::current_working_directory() << "\n";
	mFlow::WingProjectionGeometry wing;
	double aspect_ratio = 2.5; //10000
	mFlow::WingGenerators::rectangular(wing, 0.6, 3);

	mFlow::ThreeDWakeULLT sim;
	sim.wing_projection = wing;
	sim.quasi_steady = false;
	sim.vortex_ring_warping_correction = false;
	sim.symmetric = true;
	sim.symmetry_plane = HBTK::CartesianPlane(HBTK::CartesianPoint3D({ 0,0,0 }), HBTK::CartesianVector3D({ 0, 1, 0 }));
	int write_vtk_every = 10;
	bool write_inner_solutions = true;

	HBTK::AerofoilGeometry aerofoil;// = HBTK::AerofoilGenerators::sd7003();
	HBTK::CubicSpline1D camber_line = aerofoil.get_camber_spline();

	mFlow::Ramesh2014 inner_sol;
	inner_sol.camber_line = [&](double x) { return camber_line((x + 1) / 2); };
	inner_sol.camber_slope = [&](double x) {return camber_line.derivative((x + 1) / 2); };
	inner_sol.pitch_location = 1.; 
	inner_sol.delta_t = 0.00375 / 2; // dt* = 0.015
	inner_sol.free_stream_velocity.as_array() = { 0.4, 0.0 };
//	std::unique_ptr<mFlow::CanonicalFunction> heave_profile 
//		= std::make_unique<mFlow::EldredgeSmoothRamp>(0.0764 + 1, 2*0.0764 + 1, 3*0.0764 + 1, 4*0.0764 + 1, 11, 10, 0.0764, inner_sol.free_stream_velocity.magnitude());
	std::unique_ptr<mFlow::CanonicalFunction> heave_profile 
		= std::make_unique<mFlow::Harmonic>(8.0, 0.0025 * 8, 0.0);
	std::unique_ptr<mFlow::CanonicalFunction> aoa_profile_variable
		= std::make_unique<mFlow::Harmonic>(10.3, HBTK::Constants::degrees_to_radians(0.), 0.0);
	std::unique_ptr<mFlow::CanonicalFunction> aoa_const
		= std::make_unique<mFlow::ConstantValue>(HBTK::Constants::degrees_to_radians(4.)); //0.06981317);
	std::unique_ptr<mFlow::CanonicalFunction> aoa_profile
		= std::make_unique<mFlow::CanonicalFunctionSum>(std::move(aoa_const), std::move(aoa_profile_variable));
	inner_sol.foil_AoA = [&](double t)->double { return aoa_profile->f(t);  };
	inner_sol.foil_dAoAdt = [&](double t)->double { return aoa_profile->dfdx(t); };
	inner_sol.foil_Z = [&](double t)->double { return heave_profile->f(t); };
	inner_sol.foil_dZdt = [&](double t)->double { return heave_profile->dfdx(t); };
	inner_sol.number_of_fourier_terms = 8;
	inner_sol.wake_self_convection = true;

	int num_inner = 8;
	std::vector<double> inner_y_positions = HBTK::semicircspace(wing.semispan(), 0, num_inner);
	for (int i = 0; i < num_inner/2; i++) {
		sim.inner_solution_planes.push_back(HBTK::CartesianPlane(
			HBTK::CartesianPoint3D({ 0, inner_y_positions[i], 0 }),
			HBTK::CartesianPoint3D({ 1, inner_y_positions[i], 0 }),
			HBTK::CartesianPoint3D({ 0, inner_y_positions[i], 1 })));
		mFlow::Ramesh2014 inner_sol_y = inner_sol;
		inner_sol_y.semichord = wing.semichord(inner_y_positions[i]);
		sim.inner_solutions.emplace_back(inner_sol_y);
	}

	HBTK::DoubleTable table_bv, table_dwy, table_dwx, table_inner;
	table_bv.add_column("Time");
	table_bv.add_column("Step");
	table_bv.add_column("Cl");
	table_dwy.add_column("Time");
	table_dwy.add_column("Step");
	table_dwx.add_column("Time");
	table_inner.add_column("Time");
	for (int i = 0; i < (int)sim.inner_solutions.size(); i++) {
		table_bv.add_column("BV-" + std::to_string(i));
		table_dwy.add_column("DWY-" + std::to_string(i));
		table_dwx.add_column("DWX-" + std::to_string(i));
		table_inner.add_column("TE_vort-" + std::to_string(i));
	}
	for (int i = 0; i < (int)sim.inner_solutions.size(); i++) {
		table_inner.add_column("LE_vort-" + std::to_string(i));
	}

	sim.initialise();
	int n_steps = 4001;
	try {
		for (int i = 0; i < n_steps; i++) {
			std::cout << "\rStep " << i + 1 << " of " << n_steps << "        ";
			sim.advance_one_step();
			table_inner[0].push_back(i * sim.inner_solutions[0].delta_t);
			table_bv[0].push_back(i * sim.inner_solutions[0].delta_t);
			table_bv[1].push_back(i);
			table_bv[2].push_back(sim.compute_lift_coefficient());
			table_dwy[0].push_back(i * sim.inner_solutions[0].delta_t);
			table_dwy[1].push_back(i);
			table_dwx[0].push_back(i * sim.inner_solutions[0].delta_t);
			for (int ic = 0; ic < (int)sim.inner_solutions.size(); ic++) {
				table_bv[ic + 3].push_back(sim.inner_solutions[ic].bound_vorticity());
				table_inner[ic + 1].push_back(sim.inner_solutions[ic].m_te_vortex_particles.vorticity_sum());
				table_inner[ic + 1 + (int)sim.inner_solutions.size()].push_back(sim.inner_solutions[ic].m_le_vortex_particles.vorticity_sum());
				table_dwy[ic + 2].push_back(sim.inner_solutions[ic].external_purturbation(HBTK::CartesianPoint2D({ 0,0 }), sim.inner_solutions[ic].time).y());
				table_dwx[ic + 1].push_back(sim.inner_solutions[ic].external_purturbation(HBTK::CartesianPoint2D({ 0,0 }), sim.inner_solutions[ic].time).x());
			}
			if (i%write_vtk_every == 0) {
				std::ofstream wake_ostream((std::string("output/wake_") + std::to_string(i) + ".vtu").c_str());
				sim.wake_to_vtk(wake_ostream);
				if (write_inner_solutions) {
					for (int j = 0; j < (int)sim.inner_solutions.size(); j++) {
						HBTK::CartesianPlane plane = sim.inner_solution_planes[j];
						plane.origin() = plane.origin() - HBTK::CartesianVector3D({ wing.semichord(plane.origin().y()), 0, 0 });
						std::ofstream in_vort_ostream(("output/in_vort" + std::to_string(j) + "_" + std::to_string(i) + ".vtu").c_str());
						HBTK::Vtk::VtkWriter writer;
						writer.appended = false;
						writer.open_file(in_vort_ostream, HBTK::Vtk::VtkWriter::vtk_file_type::UnstructuredGrid);
						auto tev_data = sim.inner_solutions[j].m_te_vortex_particles.to_vtk_data(plane);
						auto lev_data = sim.inner_solutions[j].m_le_vortex_particles.to_vtk_data(plane);
						writer.write_piece(in_vort_ostream, tev_data);
						writer.write_piece(in_vort_ostream, lev_data);
						writer.close_file(in_vort_ostream);
					}
				}
			}
		}
	}
	catch (std::domain_error & e) {
		std::cout << "\n";
		std::cout << "Encountered error: \n";
		std::cout << e.what() << "\n";
		std::cout << "Will print any unrecorded results to file and stop.\n";
	}
	std::cout << "\n";
	HBTK::CsvWriter csv_writer;
	csv_writer.precision = 4;
	csv_writer.neat_columns = true;
	std::ofstream bv_stream("bound_vorticities.csv");
	std::ofstream dw_stream("downwash.csv");
	std::ofstream in_stream("inner.csv");
	std::ofstream fw_stream("forwash.csv");
	csv_writer.write(bv_stream, table_bv);
	csv_writer.write(dw_stream, table_dwy);
	csv_writer.write(fw_stream, table_dwx);
	csv_writer.write(in_stream, table_inner);

	return 0;
}
