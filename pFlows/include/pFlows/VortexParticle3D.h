#pragma once
/*////////////////////////////////////////////////////////////////////////////
VortexParticle3D.h

A representation a vorton / vortex particle.

Copyright 2018 HJA Bird

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

#include <cassert>
#include <cmath>

#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Constants.h>

namespace mFlow{
	class VortexParticle3D {
	public:
		VortexParticle3D() = default;
		VortexParticle3D(const HBTK::CartesianPoint3D & coordinate,
			const HBTK::CartesianVector3D & vorticity,
			const double size);
		~VortexParticle3D() = default;
		
		// Cartesian location of the vortex particle
		HBTK::CartesianPoint3D coord;

		// The vorticity vector. In the sense of the vorticity integral
		// over the particle's volume.
		HBTK::CartesianVector3D vorticity;

		// The radius of the vortex particle.
		double radius;

		HBTK::CartesianVector3D induced_velocity(const HBTK::CartesianPoint3D & measurement_point) const;

		HBTK::CartesianVector3D induced_rate_of_change_of_vorticity(
			const VortexParticle3D & other_particle
			) const;
		
	protected:
		// None!? Gee. You and your high and might OO principles,
		// but none of the encapsulation to go with it. You make
		// me sick!
	};
}


