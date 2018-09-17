#pragma once
/*////////////////////////////////////////////////////////////////////////////
VortexParticle3DCollector.h

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
