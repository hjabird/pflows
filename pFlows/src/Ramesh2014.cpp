#include "Ramesh2014.h"
/*////////////////////////////////////////////////////////////////////////////
Ramesh2014.cpp

Partial implimentation of the paper "Discrete-vortex method with novel
shedding criterion for unsteady aerofoil flows with intermittent leading-edge
vortex shedding", Kiran Ramesh et al. J. Fluid Mech 2014.
doi:10.1017/jfm.2014.297

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
#include <exception>
#include <iostream>

#include "HBTK/Checks.h"
#include "HBTK/GaussianQuadrature.h"
#include "HBTK/Generators.h"
#include "HBTK/PotentialFlowDistributions.h"


mFlow::Ramesh2014::Ramesh2014()
	: time(0),
	delta_t(0),
	semichord(0),
	pitch_location(0),
	number_of_fourier_terms(0),
	wake_self_convection(true)
{
}

mFlow::Ramesh2014::~Ramesh2014()
{
}

void mFlow::Ramesh2014::initialise()
{
	// Timestep of zero goes nowhere!
	if (delta_t <= 0) {
		throw std::domain_error("delta_t: "
			"Time step for Ramesh2014 method must be more than zero. " __FILE__
			":" + std::to_string(__LINE__));
	}

	// Timestep of zero goes nowhere!
	if (!HBTK::check_finite(delta_t)) {
		throw std::domain_error("delta_t: "
			"Time step for Ramesh2014 be finite. " __FILE__
			":" + std::to_string(__LINE__));
	}

	// But who would set time to infinite?
	if (!HBTK::check_finite(time)) {
		throw std::domain_error("time: "
			"Time value for Ramesh2014 is not finite. Cannot use non-finite"
			"values for time. " __FILE__ ":" + std::to_string(__LINE__));
	}

	// If no camber line and slope've been set, user probably wants a flat plate.
	// If not both, they've been a muppet.
	if (!camber_line || !camber_slope) {
		if (camber_line || camber_slope) {
			throw std::domain_error("Ramesh2014: "
				"Only one of the camber line and camber slope have been set - not both. "
				"If neither are set, a flat plate is assumed. Otherwise, both must be "
				"set. " __FILE__ " : " + std::to_string(__LINE__));
		}
		else {
			camber_line = [](double x)->double { return 0; };
			camber_slope = [](double x)->double { return 0; };
		}
	}

	// Fourier terms:
	if (number_of_fourier_terms == 0) {
		throw std::domain_error("number_of_fourier_terms: "
			"number of fourier terms for Ramesh2014 method should be more than 2. " 
			__FILE__ ":"  + std::to_string(__LINE__));
	}
	m_fourier_terms = HBTK::uniform(0.0, number_of_fourier_terms);
	time -= delta_t;
	compute_fourier_terms();
	time += delta_t;
	m_previous_fourier_terms = m_fourier_terms;
}

void mFlow::Ramesh2014::advance_one_step()
{
	convect_particles();
	shed_new_particle();
	adjust_last_shed_vortex_particle_for_kelvin_condition();
	m_previous_fourier_terms = m_fourier_terms;
	compute_fourier_terms();
	calculate_velocities();
	time += delta_t;
}

int mFlow::Ramesh2014::number_of_particles()
{
	return (int) m_vortex_particles.size();
}


std::vector<double> mFlow::Ramesh2014::get_fourier_terms()
{
	return m_fourier_terms;
}

void mFlow::Ramesh2014::calculate_velocities()
{
	// Initialise at zero for sum.
	for (auto & particle_vel : m_vortex_particle_velocities) {
		particle_vel.x() = 0;
		particle_vel.y() = 0;
	}

	// Velocities induced on particles by other particles (N-Body problem)
	if(wake_self_convection){
		double vc_size = vortex_core_size();
		for (int i = 0; i < number_of_particles(); i++) {
			double xi, yi;	// Cartesian coordinates of ith particle.
			xi = m_vortex_particles[i].position.x();
			yi = m_vortex_particles[i].position.y();
			for (int j = 0; j < i; j++) {
				double u, v;	// Cartesian velocity components
				double xj, yj;	// Cartesian coordinates of jth particle
				xj = m_vortex_particles[j].position.x();
				yj = m_vortex_particles[j].position.y();
				std::tie(u, v) = unity_vortex_blob_induced_vel(xj, yj, xi, yi, vc_size);
				m_vortex_particle_velocities[j].x() += u * m_vortex_particles[i].vorticity;
				m_vortex_particle_velocities[j].y() += v * m_vortex_particles[i].vorticity;
				m_vortex_particle_velocities[i].x() -= u * m_vortex_particles[j].vorticity;
				m_vortex_particle_velocities[i].y() -= v * m_vortex_particles[j].vorticity;
			}
		}
	}

	// Velocity induced on particles by the aerofoil's vorticity distribution.
	{
		// We'll integrate over the aerofoil's vorticity distribution using a remapped Gauss quadrature.
		HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
		quad.telles_cubic_remap(-1.0);	// Adapt quadrature for the leading edge singularity.
		auto u_integrand_outer = [&](double local_foil_pos, double x_mes, double y_mes)->double {
			double xf, yf;	// Cartesian coordinates of the corresponing point on the aerofoil.
			std::tie(xf, yf) = foil_coordinate(local_foil_pos);
			double u = HBTK::PointVortex::unity_u_vel(x_mes, y_mes, xf, yf) 
				* vorticity_density(local_foil_pos) * semichord;
			return u;
		};
		auto v_integrand_outer = [&](double local_foil_pos, double x_mes, double y_mes)->double {
			double xf, yf;	// Cartesian coordinates of the corresponing point on the aerofoil.
			std::tie(xf, yf) = foil_coordinate(local_foil_pos);
			double v = HBTK::PointVortex::unity_v_vel(x_mes, y_mes, xf, yf) 
				* vorticity_density(local_foil_pos) * semichord;
			return v;
		};

		for (int i = 0; i < number_of_particles(); i++) {
			double x_mes = m_vortex_particles[i].position.x();
			double y_mes = m_vortex_particles[i].position.y();
			auto u_integrand_inner = [&](double local_pos) { return u_integrand_outer(local_pos, x_mes, y_mes); };
			auto v_integrand_inner = [&](double local_pos) { return v_integrand_outer(local_pos, x_mes, y_mes); };
			m_vortex_particle_velocities[i].x() += quad.integrate(u_integrand_inner);
			m_vortex_particle_velocities[i].y() += quad.integrate(v_integrand_inner);
		}
	}

	// Velocity induced by free stream.
	for (auto & particle_vel : m_vortex_particle_velocities) {
		particle_vel.x() += free_stream_velocity.x();
		particle_vel.y() += free_stream_velocity.y();
	}

	// We're done. Check that we've done something sane:
	for (auto particle_vel : m_vortex_particle_velocities) {
		assert(HBTK::check_finite(particle_vel.x()));
		assert(HBTK::check_finite(particle_vel.y()));
	}
}

void mFlow::Ramesh2014::convect_particles()
{
	// Forward Euler scheme - justified by Ramesh et al.
	for (int i = 0; i < number_of_particles(); i++) {
		m_vortex_particles[i].position.x() += m_vortex_particle_velocities[i].x() * delta_t;
		m_vortex_particles[i].position.y() += m_vortex_particle_velocities[i].y() * delta_t;
	}
	return;
}

void mFlow::Ramesh2014::shed_new_particle()
{
	VortexGroup2D::Vortex2D particle;
	HBTK::CartesianVector2D velocity;
	if (number_of_particles() > 0) {
		particle = m_vortex_particles.most_recently_added();
		double te_x, te_y;	// Cartesian coordinates of trailing edge.
		std::tie(te_x, te_y) = foil_coordinate(1);
		// Adjust so that the new particle is 1/3 of the distance
		// between the TE and the last shed particle.
		particle.position.x() -= (2. / 3.) * (particle.position.x() - te_x);
		particle.position.y() -= (2. / 3.) * (particle.position.y() - te_y);
		particle.vorticity = 0;
	}
	else {
		// Our first particle.
		// Put the particle in the path of foil.
		particle.vorticity = 0;
		std::tie(particle.position.x(), particle.position.y()) = foil_coordinate(1);
		std::tie(velocity.x(), velocity.y()) = foil_velocity(1);
		velocity = velocity - free_stream_velocity;
		// 0.5 corresponds to the 1/3 rule used earlier (0.5 + 1 = 3 * 0.5)
		particle.position.x() -= velocity.x() * delta_t * 0.5;
		particle.position.y() -= velocity.y() * delta_t * 0.5;
	}
	m_vortex_particles.add(particle);
	m_vortex_particle_velocities.emplace_back(velocity);
	return;
}

double mFlow::Ramesh2014::total_shed_vorticity()
{
	double vorticity = 0;
	for (int i = 0; i < m_vortex_particles.size(); i++) {
		vorticity += m_vortex_particles[i].vorticity;
	}
	if (!HBTK::check_finite(vorticity)) {
		std::overflow_error("Infinite vorticity: "
			"The vorticity of the wake was found to be infinite whilst "
			"running Ramesh2014. " __FILE__ ":" + std::to_string(__LINE__));
	}
	return vorticity; 
}

std::pair<double, double> mFlow::Ramesh2014::get_particle_induced_velocity(
	double x, double y)
{
	double vx = 0;	// Cartesian x velocity
	double vy = 0;	// Cartesian y velocity
	double vc_size = vortex_core_size();
	for (int i = 0; i < m_vortex_particles.size(); i++) {
		double u, v;	// Velocity induced by single vortex blob.
		std::tie(u, v) = unity_vortex_blob_induced_vel(x, y,
			m_vortex_particles[i].position.x(), m_vortex_particles[i].position.y(), vc_size);
		vx += m_vortex_particles[i].vorticity * u;
		vy += m_vortex_particles[i].vorticity * v;
		assert(HBTK::check_finite(vx));
		assert(HBTK::check_finite(vy));
	}
	return std::make_pair(vx, vy);
}


void mFlow::Ramesh2014::adjust_last_shed_vortex_particle_for_kelvin_condition()
{
	m_vortex_particles.most_recently_added().vorticity = 0;
	// The bound vorticity distribution can be decomposed into a part that is a result of
	// the wake for which we have already set the vorticity.
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	auto known_integrand = [&](double theta) -> double {
		double local_x = - cos(theta);		// [0, pi] -> x in [-1,1] -> [LE, TE]
		double wake_induced_v, wake_induced_u;
		double x, y;	// Coordinate of point we're evaluating
		std::tie(x, y) = foil_coordinate(local_x);
		std::tie(wake_induced_u, wake_induced_v) = get_particle_induced_velocity(x, y);
		double local_camber_slope = camber_slope(local_x);
		double T_1 =  local_camber_slope *
			(free_stream_velocity.x() * cos(alpha) 
				+ free_stream_velocity.y() * sin(alpha)
				+ h_dot * sin(alpha)
				- sin(alpha) * wake_induced_u  
				+ cos(alpha) * wake_induced_u )
			- free_stream_velocity.x() * sin(alpha)
			- free_stream_velocity.y() * cos(alpha)
			- alpha_dot * semichord * (1 + local_x - 2 * pitch_location)
			+ h_dot * cos(alpha)
			- cos(alpha) * wake_induced_v
			- sin(alpha) * wake_induced_u;
		return T_1 * (cos(theta) - 1) * 2.0 * semichord;
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	double I_known = quad.integrate(known_integrand);
	assert(HBTK::check_finite(I_known));

	// And a part that will be cause by our new vortex particle.
	double I_unknown;
	{
		double x_p, y_p; // Position of our last shed vortex particle.
		x_p = m_vortex_particles.most_recently_added().position.x();
		y_p = m_vortex_particles.most_recently_added().position.y();
		double v_core_size = vortex_core_size();
		auto particle_integrand = [&](double theta) -> double {
			double foil_pos = - cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
			double x_f, y_f; // Position on our foil.
			std::tie(x_f, y_f) = foil_coordinate(foil_pos);
			double vx, vy;
			std::tie(vx, vy) = unity_vortex_blob_induced_vel(x_f, y_f, x_p, y_p, v_core_size);
			double normal_velocity = -cos(alpha) * vy - sin(alpha) * vx;
			return normal_velocity * (cos(theta) - 1) * 2.0 * semichord;
		};
		// We expect particle_integrand to be singular looking towards the TE.
		quad.telles_quadratic_remap(0);	
		I_unknown = quad.integrate(particle_integrand);
	}
	m_vortex_particles.most_recently_added().vorticity = - (I_known + total_shed_vorticity()) / (1 + I_unknown);
	return;
}

void mFlow::Ramesh2014::compute_fourier_terms()
{
	assert(number_of_fourier_terms > 2);
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	for (int i = 0; i < number_of_fourier_terms; i++) {
		auto integrand = [&](double theta) {
			double foil_pos = - cos(theta);
			double x, y;		// Coordinate (global)
			std::tie(x, y) = foil_coordinate(foil_pos);
			double pvx, pvy;	// vortex particle induced velocities
			std::tie(pvx, pvy) = get_particle_induced_velocity(x, y);
			double alpha = foil_AoA(time);
			double alpha_dot = foil_dAoAdt(time);
			double h_dot = foil_dZdt(time);
			// Return velocity normal to the foil.
			double local_camber_slope = camber_slope(foil_pos);
			return cos(i * theta) * (local_camber_slope *
				(free_stream_velocity.x() * cos(alpha)
					+ free_stream_velocity.y() * sin(alpha)
					+ h_dot * sin(alpha)
					- sin(alpha) * pvy 
					+ cos(alpha) * pvx)
				- free_stream_velocity.x() * sin(alpha)
				- free_stream_velocity.y() * cos(alpha)
				- alpha_dot * semichord * (1 - cos(theta) - 2 * pitch_location)
				+ h_dot * cos(alpha)
				- cos(alpha) * pvy - sin(alpha) * pvx)
				/ free_stream_velocity.magnitude(); 
		};
		m_fourier_terms[i] = quad.integrate(integrand) * 2 / HBTK::Constants::pi();
	}
	m_fourier_terms[0] *= -0.5;
	return;
}

double mFlow::Ramesh2014::bound_vorticity()
{	// Eq.2.7
	return free_stream_velocity.magnitude() *
		semichord * HBTK::Constants::pi() * 
		(2 * m_fourier_terms[0] + m_fourier_terms[1]);
}

double mFlow::Ramesh2014::vorticity_density(double local_pos)
{
	assert(HBTK::check_finite(m_fourier_terms));
	assert(local_pos != -1.);	// Leading edge singularity will cause problems.
	double theta = acos(-local_pos);
	double vd = (local_pos == 1. ? 0 : m_fourier_terms[0] * (1 + cos(theta)) / sin(theta));
	for (int i = 1; i < (int)m_fourier_terms.size(); i++) {
		vd += m_fourier_terms[i] * sin(i * theta);
	}
	vd *= 2 * free_stream_velocity.magnitude();
	return vd;
}

std::pair<double, double> mFlow::Ramesh2014::foil_coordinate(double eta)
{
	double x, y;
	double alpha = foil_AoA(time);
	x = semichord * eta;
	y = camber_line(eta);
	x -= semichord * pitch_location;
	double xtmp = x * cos(alpha) + y * sin(alpha);
	double ytmp = -x * sin(alpha) + y * cos(alpha);
	x = xtmp + semichord * pitch_location;
	y = ytmp + foil_Z(time);
	return std::make_pair(x, y);
}

std::pair<double, double> mFlow::Ramesh2014::foil_velocity(double eta)
{
	double angular_vel = foil_dAoAdt(time);
	double AoA = foil_AoA(time);
	double x, y, xp, yp;
	std::tie(x, y) = foil_coordinate(eta);
	std::tie(xp, yp) = foil_coordinate(0);
	double radius = sqrt(pow(x-xp, 2) + pow(y-yp, 2));
	double angle = atan2(yp, xp);
	double vy = angular_vel * radius * cos(AoA);
	double vx = angular_vel * radius * sin(AoA);
	vy += foil_dZdt(time);
	return std::make_pair(vx, vy);
}

double mFlow::Ramesh2014::aerofoil_leading_edge_suction_force(double density)
{
	assert(HBTK::check_finite(m_fourier_terms));
	return pow(free_stream_velocity.x(), 2) * 
		pow(free_stream_velocity.y(), 2) * HBTK::Constants::pi() *
		density * 2 * semichord * pow(m_fourier_terms[0], 2);
}

double mFlow::Ramesh2014::aerofoil_normal_force(double density)
{
	assert(HBTK::check_finite(m_fourier_terms));
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	std::vector<double> fourier_derivatives = rate_of_change_of_fourier_terms();

	double term_1, term_2,
		term_11, term_12, term_121, term_122,
		term_1211, term_1212;

	term_11 = density * HBTK::Constants::pi() * semichord * 2.0 * free_stream_velocity.magnitude();
	term_1211 = free_stream_velocity.x() * cos(alpha) 
		+ free_stream_velocity.y() * sin(alpha) + h_dot * sin(alpha);
	term_1212 = m_fourier_terms[0] + m_fourier_terms[1] / 2;
	term_121  = term_1211 * term_1212;
	term_122 = 2 * semichord * (
		(3. / 4) * fourier_derivatives[0]
		+ (1. / 4) * fourier_derivatives[1]
		+ (1. / 8) * fourier_derivatives[2]);
	term_12 = term_121 + term_122;
	term_1 = term_11 * term_12;
	assert(HBTK::check_finite(term_1));

	// term_2 includes an integral.
	auto integrand = [&](double local_pos)->double {
		double vort_den = vorticity_density(local_pos);
		double x, y;
		std::tie(x, y) = foil_coordinate(local_pos);
		double u, w; 
		std::tie(u, w) = get_particle_induced_velocity(x, y);
		double tangential_velocity = u * cos(alpha) - w * sin(alpha);
		return tangential_velocity * vort_den;
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	term_2 = quad.integrate(integrand) * semichord * density;
	assert(HBTK::check_finite(term_2));
	return term_1 + term_2;
}

double mFlow::Ramesh2014::aerofoil_moment_about_pitch_location(double density)
{
	assert(HBTK::check_finite(density));
	assert(HBTK::check_finite(m_fourier_terms));
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	std::vector<double> fourier_derivatives = rate_of_change_of_fourier_terms();

	double term_1, term_2, term_3,
		term_21, term_22, term_211, term_212;
	term_1 = aerofoil_normal_force(density) * (1 + pitch_location) * semichord;
	term_211 = free_stream_velocity.x() * cos(alpha)
		+ free_stream_velocity.y() * sin(alpha) + h_dot * sin(alpha);
	term_212 = 0.25 * m_fourier_terms[0] + 0.25 * m_fourier_terms[1]
		- 0.25 * m_fourier_terms[2];
	term_21 = term_211 * term_212;
	term_22 = semichord * 2 * (
		(7. / 16) * fourier_derivatives[0]
		+ (11. / 64) * fourier_derivatives[1]
		+ (1. / 16) * fourier_derivatives[2]
		- (1. / 64) * fourier_derivatives[3]);
	term_2 = (term_21 + term_22) * density * 4 * 
		free_stream_velocity.magnitude() * pow(semichord, 2);
	assert(HBTK::check_finite(term_2));

	// term_3 includes an integral.
	auto integrand = [&](double local_pos)->double {
		double vort_den = vorticity_density(local_pos);
		double x, y;
		std::tie(x, y) = foil_coordinate(local_pos);
		double u, w;
		std::tie(u, w) = get_particle_induced_velocity(x, y);
		double tangential_velocity = u * cos(alpha) - w * sin(alpha);
		return tangential_velocity * vort_den * semichord * (1 + local_pos);
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(40);
	term_3 = quad.integrate(integrand) * semichord * density;
	assert(HBTK::check_finite(term_3));

	return term_1 - term_2 - term_3;
}

std::pair<double, double> mFlow::Ramesh2014::aerofoil_lift_and_drag_coefficients()
{
	double stream_vel_norm = pow(free_stream_velocity.x(), 2);
	double normal_force_coeff = aerofoil_normal_force(1.) / 
		(stream_vel_norm * semichord);
	double suction_force_coeff = aerofoil_leading_edge_suction_force(1.) /
		(stream_vel_norm * semichord);
	double cl, cd;
	double alpha = foil_AoA(time);
	cl = normal_force_coeff * cos(alpha) + suction_force_coeff * sin(alpha);
	cd = normal_force_coeff * sin(alpha) - suction_force_coeff * cos(alpha);
	return std::make_pair(cl, cd);
}

std::vector<double> mFlow::Ramesh2014::rate_of_change_of_fourier_terms()
{
	assert(delta_t != 0);
	std::vector<double> time_derivatives(m_fourier_terms.size());
	for (int i = 0; i < (int) m_fourier_terms.size(); i++) {
		time_derivatives[i] = (m_fourier_terms[i] - m_previous_fourier_terms[i]) / delta_t;
	}
	assert(HBTK::check_finite(time_derivatives));
	return time_derivatives;
}

double mFlow::Ramesh2014::vortex_core_size() const
{
	double delta_t_star_times_c = delta_t * free_stream_velocity.magnitude();
	return delta_t_star_times_c * 1.3;
}

std::pair<double, double> mFlow::Ramesh2014::unity_vortex_blob_induced_vel(
	double x_mes, double y_mes, double x_vor, double y_vor, double vortex_size)
{
	assert((x_mes != x_vor) || (y_mes != y_vor));
	double u, v, denominator;
	denominator = sqrt(pow(pow(x_mes - x_vor, 2) + pow(y_mes - y_vor, 2), 2)
		+ pow(vortex_size, 4)) * 2. * HBTK::Constants::pi();
	u = (y_mes - y_vor) / denominator;
	v = - (x_mes - x_vor) / denominator;
	return std::make_pair(u, v);
}
