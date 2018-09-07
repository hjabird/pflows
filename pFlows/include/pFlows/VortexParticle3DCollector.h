#pragma once

#include <fstream>
#include <vector>

#include "AbstractVortexObject.h"
#include "VortexParticle3D.h"
#include "VortexParticleKernels.h"

namespace mFlow {
	class VortexParticle3DCollector
		: AbstractVortexObject
	{
	public:
		std::vector<VortexParticle3D> particles;
		VortexParticleRegularisation regularisation;

		~VortexParticle3DCollector();

		// The induced velocity at any point:
		virtual HBTK::CartesianVector3D induced_velocity(const HBTK::CartesianPoint3D& mes_pnt);

		// The rate of change of vorticity on a vortex particle.
		virtual HBTK::CartesianVector3D induced_dvort(const VortexParticle3D& particle);

		// Add the geometry and the vorticity to a vtk unstructured data set.
		virtual void add_to_vtk(HBTK::Vtk::VtkUnstructuredDataset& dataset);

	protected:
	};
}
