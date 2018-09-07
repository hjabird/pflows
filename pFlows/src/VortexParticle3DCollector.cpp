#include "VortexParticle3DCollector.h"

#include <HBTK/VtkUnstructuredMeshHolder.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

mFlow::VortexParticle3DCollector::~VortexParticle3DCollector()
{
}

HBTK::CartesianVector3D mFlow::VortexParticle3DCollector::induced_velocity(const HBTK::CartesianPoint3D & mes_pnt)
{
	double u(0), v(0), w(0);
#pragma omp parallel for reduction(+ : u), reduction(+ : v), reduction(+ : w)
	for (int i = 0; i < (int)particles.size(); i++) {
		VortexParticle3D& p = particles[i];
		HBTK::CartesianVector3D vel = p.induced_velocity(mes_pnt);
		u += vel.x();
		v += vel.y();
		w += vel.z();
	}
	return HBTK::CartesianVector3D({ u, v, w });
}

HBTK::CartesianVector3D mFlow::VortexParticle3DCollector::induced_dvort(const VortexParticle3D & particle)
{
	double u(0), v(0), w(0);
#pragma omp parallel for reduction(+ : u), reduction(+ : v), reduction(+ : w)
	for (int i = 0; i < (int)particles.size(); i++) {
		VortexParticle3D& p = particles[i];
		HBTK::CartesianVector3D vel = p.induced_rate_of_change_of_vorticity(particle);
		u += vel.x();
		v += vel.y();
		w += vel.z();
	}
	return HBTK::CartesianVector3D({ u, v, w });
}

void mFlow::VortexParticle3DCollector::add_to_vtk(HBTK::Vtk::VtkUnstructuredDataset & dataset)
{
	HBTK::Vtk::VtkUnstructuredDataset & data = dataset;
	HBTK::Vtk::VtkUnstructuredMeshHolder & mesh = data.mesh;

	for (int i = 0; i < (int)particles.size(); i++) {
		mesh.points.push_back(particles[i].coord);
		mesh.cells.push_back({ 1, std::vector<int>({ i }) });
		data.vector_cell_data["vorticity"].push_back(particles[i].vorticity);
	}
	return;
}

