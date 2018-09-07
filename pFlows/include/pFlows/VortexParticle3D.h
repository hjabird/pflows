#pragma once

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
		
		HBTK::CartesianPoint3D coord;
		HBTK::CartesianVector3D vorticity;
		double size;

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


