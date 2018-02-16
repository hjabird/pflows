#include "stdafx.h"
#include "Ramesh2014.h"

#include <cassert>
#include "HBTK/PotentialFlowDistributions.h"
#include "HBTK/Checks.h"
#include "HBTK/GaussianQuadrature.h"


mFlow::Ramesh2014::Ramesh2014()
	: time(0),
	delta_t(0),
	semichord(0),
	pitch_location(0)
{
}

mFlow::Ramesh2014::~Ramesh2014()
{
}

int mFlow::Ramesh2014::num_particles()
{
	return (int) m_vortex_particles.size();
}


void mFlow::Ramesh2014::calculate_velocities()
{
	// Initialise at zero for sum.
	for (auto & particle : m_vortex_particles) {
		particle.vx = 0;
		particle.vy = 0;
	}

	// particle velocities.
	for (int i = 0; i < num_particles(); i++) {
		double xi, yi;
		xi = m_vortex_particles[i].x;
		yi = m_vortex_particles[i].y;
		for (int j = 0; j < i; j++) {
			double u, v;
			double xj, yj;
			xj = m_vortex_particles[j].x;
			yj = m_vortex_particles[j].y;
			u = HBTK::PointVortex::unity_u_vel(xj, yj, xi, yi);
			v = HBTK::PointVortex::unity_v_vel(xj, yj, xi, yi);
			assert(HBTK::check_finite(u));
			assert(HBTK::check_finite(v));
			m_vortex_particles[j].vx += u * m_vortex_particles[i].vorticity;
			m_vortex_particles[j].vy += v * m_vortex_particles[i].vorticity;
			m_vortex_particles[i].vx -= u * m_vortex_particles[j].vorticity;
			m_vortex_particles[i].vy -= v * m_vortex_particles[j].vorticity;
		}
	}

	// Foil induced velocities.


}

void mFlow::Ramesh2014::convect_particles()
{
	for (auto & particle : m_vortex_particles) {
		particle.x += particle.vx * delta_t;
		particle.y += particle.vy * delta_t;
	}
	return;
}

void mFlow::Ramesh2014::shed_new_particle()
{
	vortex_particle particle;
	if (num_particles() > 0) {
		particle.vorticity = 0;
		std::tie(particle.x, particle.y) = foil_coordinate(1);
		std::tie(particle.vx, particle.vy) = foil_velocity(1);
		particle.vx += free_stream_velocity;
		particle.x += particle.vx * delta_t * 0.5;
		particle.y += particle.vy * delta_t * 0.5;
	}
	else {
		particle = m_vortex_particles.back();
		double te_x, te_y;
		std::tie(te_x, te_y) = foil_coordinate(1);
		particle.x -= (2. / 3.) * (particle.x - te_x);
		particle.y -= (2. / 3.) * (particle.y - te_y);
	}
	m_vortex_particles.emplace_back(particle);
	return;
}

double mFlow::Ramesh2014::shed_vorticity()
{
	double vorticity = 0;
	for (auto particle : m_vortex_particles) {
		vorticity += particle.vorticity;
	}
	assert(HBTK::check_finite(vorticity));
	return vorticity; 
}

std::pair<double, double> mFlow::Ramesh2014::get_particle_induced_velocity(double x, double y)
{
	double vx = 0;
	double vy = 0;
	for (auto particle : m_vortex_particles) {
		vx += particle.vorticity *
			HBTK::PointVortex::unity_u_vel(x, y, particle.x, particle.y);
		vy += particle.vorticity *
			HBTK::PointVortex::unity_v_vel(x, y, particle.x, particle.y);
		assert(HBTK::check_finite(vx));
		assert(HBTK::check_finite(vy));
	}
	return std::make_pair(vx, vy);
}


void mFlow::Ramesh2014::adjust_last_shed_vortex_particle_for_kelvin_condition()
{
	double system_vorticity = shed_vorticity() + bound_vorticity();
	assert(num_particles() > 0);
	m_vortex_particles.back().vorticity = -system_vorticity;
	return;
}

void mFlow::Ramesh2014::compute_fourier_terms()
{
	assert(num_fourier_terms > 2);
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(40);
	quad.linear_remap(0, HBTK::Constants::pi());
	for (int i = 0; i < num_fourier_terms; i++) {
		auto integrand = [&](double theta) {
			double foil_pos = semichord * (1 - cos(theta));
			double x, y;		// Coordinate (global)
			std::tie(x, y) = foil_coordinate(foil_pos);
			double pvx, pvy;	// vortex particle induced velocities
			std::tie(pvx, pvy) = get_particle_induced_velocity(x, y);
			double alpha = foil_AoA(time);
			double alpha_dot = foil_dAoAdt(time);
			double h_dot = foil_dZdt(time);
			double fvx, fvy;	// Foil surface velocity.
			std::tie(fvx, fvy) = foil_velocity(foil_pos);
			double rvx, rvy;    // Relative velocity of foil vortex particle induced flow.
			rvx = pvx - fvx;
			rvy = pvy - fvy;
			// Return velocity normal to the foil.
			return cos(i * theta) * (- rvx * sin(alpha) + rvy * cos(alpha)) 
				/ free_stream_velocity;
		};
		m_fourier_terms[i] = quad.integrate(integrand) * 2 / HBTK::Constants::pi();
	}
	m_fourier_terms[0] *= -0.5;
	return;
}

double mFlow::Ramesh2014::bound_vorticity()
{
	return free_stream_velocity * semichord * HBTK::Constants::pi() * 
		(2 * m_fourier_terms[0] + m_fourier_terms[1]);
}

std::pair<double, double> mFlow::Ramesh2014::foil_coordinate(double eta)
{
	double x, y;
	x = semichord * (pitch_location - eta) * cos(foil_AoA(time));
	y = foil_Z(time) + semichord * (pitch_location - eta) * sin(foil_AoA(time));
	return std::make_pair(x, y);
}

std::pair<double, double> mFlow::Ramesh2014::foil_velocity(double eta)
{
	double vx, vy;
	vy = foil_dZdt(time);
	double radial_vel = semichord * (pitch_location - eta) * foil_dAoAdt(time);
	vy += radial_vel * sin(foil_AoA(time));
	vx = radial_vel * cos(foil_AoA(time));
	return std::make_pair(vx, vy);
}
