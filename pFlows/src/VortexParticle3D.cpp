#include "VortexParticle3D.h"

#include <cassert>

#include <HBTK/Checks.h>
#include <HBTK/CartesianVector.h>

#include "VortexParticleKernels.h"

mFlow::VortexParticle3D::VortexParticle3D(
	const HBTK::CartesianPoint3D & coordinate, 
	const HBTK::CartesianVector3D & vorticity, 
	const double size)
	:	coord(coordinate),
	vorticity(vorticity),
	radius(size)
{
	assert(HBTK::check_finite(coordinate));
	assert(HBTK::check_finite(vorticity));
	assert(HBTK::check_finite(size));
	assert(size > 0);
}

HBTK::CartesianVector3D mFlow::VortexParticle3D::induced_velocity(
	const HBTK::CartesianPoint3D & measurement_point) const
{
	assert(HBTK::check_finite(measurement_point));
	HBTK::CartesianVector3D vel({ 0,0,0 });
	if (measurement_point != coord) {
		double term1 = 1. / (4. * HBTK::Constants::pi());
		HBTK::CartesianVector3D rad = coord - measurement_point;
		double term2d = pow(abs(rad), 3);
		double rho = radius / abs(rad);
		double g = winckelmans_vortex::g_function(rho);
		HBTK::CartesianVector3D term2n = g * rad.cross(vorticity);
		HBTK::CartesianVector3D term2 = term2n / term2d;
		vel = term1 * term2;
	}
	return vel;
}


HBTK::CartesianVector3D mFlow::VortexParticle3D::induced_rate_of_change_of_vorticity(
	const VortexParticle3D & other_particle) const
{
	HBTK::CartesianVector3D dvort({ 0,0,0 });

	if (&other_particle != this) {
		HBTK::CartesianVector3D r = coord - other_particle.coord;
		double sigma_k = other_particle.radius;
		double rho = abs(r) / radius;
		double g = winckelmans_vortex::g_function(rho);
		double f = winckelmans_vortex::f_function(rho);

		double term1 = 1. / (4. * HBTK::Constants::pi() * pow(sigma_k, 3));
		HBTK::CartesianVector3D  term21 = -g * vorticity.cross(other_particle.vorticity) / pow(abs(r) / sigma_k, 3);
		double term221 = 1. / pow(abs(r), 2);
		double term222 = 3 * g / pow(abs(r) / sigma_k, 3) - f;
		HBTK::CartesianVector3D term223 = vorticity.dot(r) * r.cross(other_particle.vorticity);

		dvort = term1 * (term21 + term221 * term222 * term223);
	}
	return dvort;
}