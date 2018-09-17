#include "VortexParticle3DCollector.h"
/*////////////////////////////////////////////////////////////////////////////
VortexParticle3DCollector.cpp

An object to hold a collection of vortex particles. Derived from
AbstractVortexObject.

Copyright 2017-2018 HJA Bird

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

