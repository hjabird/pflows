#pragma once

#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/VtkUnstructuredDataset.h>

#include "VortexParticle3D.h"

namespace mFlow {
	class AbstractVortexObject {
	public:
		// Construct the object
		//virtual ~AbstractVortexObject() = 0;

		// The induced velocity at any point:
		virtual HBTK::CartesianVector3D induced_velocity(const HBTK::CartesianPoint3D& mes_pnt) = 0;

		// The rate of change of vorticity on a vortex particle.
		virtual HBTK::CartesianVector3D induced_dvort(const VortexParticle3D& particle) = 0;

		// Add the geometry and the vorticity to a vtk unstructured data set.
		virtual void add_to_vtk(HBTK::Vtk::VtkUnstructuredDataset& dataset) = 0;
	};
}