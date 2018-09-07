
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>

#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianPlane.h>
#include <HBTK/Constants.h>
#include <HBTK/Paths.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

#include "VortexParticle3D.h"
#include "VortexParticle3DCollector.h"



int main(int argc, char* argv[]) {
	std::cout << "Vortex particle demo (C) HJA Bird 2018\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";

	double dt = 0.03;
	int n_steps = 2000;
	int write_every = 200;

	int n_per_ring = 30;
	double r = 1;
	double vorticity = 1;
	HBTK::CartesianPlane plane1(
		HBTK::CartesianPoint3D({ 0., 0., 0. }),
		HBTK::CartesianPoint3D({ 1., 0., 0. }),
		HBTK::CartesianPoint3D({ 0., 1., 0. })
	);
	HBTK::CartesianPlane plane2(
		HBTK::CartesianPoint3D({ 0., 0., 1. }),
		HBTK::CartesianPoint3D({ 1., 0., 1. }),
		HBTK::CartesianPoint3D({ 0., 1., 1. })
	);

	mFlow::VortexParticle3DCollector particles;
	// Make a ring.
	for (int i = 0; i < n_per_ring; ++i) {
		double dtheta = 2 * HBTK::Constants::pi() / n_per_ring;
		double theta = i * dtheta;
		HBTK::CartesianPoint2D pnt({ r * cos(theta), r * sin(theta) });
		mFlow::VortexParticle3D p;
		p.coord = plane1(pnt);
		HBTK::CartesianVector3D a, b, c;
		a = plane1(pnt) - plane1.origin();
		b = plane1.normal();
		c = a.cross(b);
		c.normalise();
		c *= dtheta * r * vorticity;
		p.vorticity = c;
		p.size = dtheta * r;
		particles.particles.push_back(p);
	}
	for (int i = 0; i < n_per_ring; ++i) {
		double dtheta = 2 * HBTK::Constants::pi() / n_per_ring;
		double theta = i * dtheta;
		HBTK::CartesianPoint2D pnt({ r * cos(theta), r * sin(theta) });
		mFlow::VortexParticle3D p;
		p.coord = plane2(pnt);
		HBTK::CartesianVector3D a, b, c;
		a = plane2(pnt) - plane2.origin();
		b = plane2.normal();
		c = a.cross(b);
		c.normalise();
		c *= dtheta * r * vorticity;
		p.vorticity = c;
		p.size = dtheta * r;
		particles.particles.push_back(p);
	}

	std::cout << "Integrating in time...\n";
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n_steps; ++i){
		std::vector<HBTK::CartesianVector3D> vels, dvorts;
		for (auto& p : particles.particles) {
			vels.push_back(particles.induced_velocity(p.coord));
			dvorts.push_back(particles.induced_dvort(p));
		}
		for (int i = 0; i < (int)particles.particles.size(); i++) {
			auto & p = particles.particles[i];
			p.coord += dt * vels[i];
			p.vorticity += dt * dvorts[i];
		}

		if (i % write_every == 0) {
			HBTK::Vtk::VtkUnstructuredDataset data;
			particles.add_to_vtk(data);
			HBTK::Vtk::VtkWriter writer;
			writer.ascii = false;
			writer.appended = false;
			std::ofstream file(std::string("output/vtk_data_" + std::to_string(i) + ".vtu").c_str());
			writer.open_file(file, HBTK::Vtk::VtkWriter::vtk_file_type::UnstructuredGrid);
			writer.write_piece(file, data);
			writer.close_file(file);
		}
	}
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "elapsed time: " << time.count() << " seconds\n";
	return 0;
}
