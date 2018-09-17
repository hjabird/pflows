#pragma once
/*////////////////////////////////////////////////////////////////////////////
AbstractVortexObject.h

An abstract base class for vortex objects: objects that represent vorticity
in boundary element sense.

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