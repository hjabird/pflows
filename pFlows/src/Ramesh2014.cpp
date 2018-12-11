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
#include "HBTK/Integrators.h"
#include "HBTK/PotentialFlowDistributions.h"


mFlow::Ramesh2014::Ramesh2014()
	: time(0),
	delta_t(0),
	semichord(0),
	pitch_location(0),
	number_of_fourier_terms(0),
	wake_self_convection(true),
	lev_shedding(false),
	shed_zero_strength_levs(false),
	critical_leading_edge_suction(1e99),
	m_wake_vorticity(0),
	m_shed_lev_at_last_step(false),
	external_purturbation([](HBTK::CartesianPoint2D p, double t) { return HBTK::CartesianVector2D({ 0,0 }); })
{
	std::cout << "FYI, I've fiddled with advance_one_step so it isn't valid for LDVM right now.\n";
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
	if (number_of_fourier_terms < 4) {
		throw std::domain_error("number_of_fourier_terms: "
			"number of fourier terms for Ramesh2014 method should be more than 3. " 
			__FILE__ ":"  + std::to_string(__LINE__));
	} 
	m_fourier_terms = HBTK::uniform(0.0, number_of_fourier_terms);
	time -= delta_t;
	compute_fourier_terms();
	time += delta_t;
	m_previous_fourier_terms = m_fourier_terms;
	m_wake_vorticity = total_shed_vorticity();
}

void mFlow::Ramesh2014::advance_one_step()
{
	convect_particles();
	time += delta_t;
	m_downwash_cache.clear();
	// m_previous_fourier_terms = m_fourier_terms; // "Wrong" location for Ramesh2014, better results for 15...
	shed_new_trailing_edge_particle_with_zero_vorticity();
	adjust_last_shed_vortex_particle_for_kelvin_condition();
	compute_fourier_terms();
	calc_rate_of_change_of_fourier_terms();
	shed_new_leading_edge_particle_if_required_and_adjust_vorticities();
	compute_fourier_terms();
	m_previous_fourier_terms = m_fourier_terms;
	calculate_velocities();
	m_downwash_cache.clear();
}

int mFlow::Ramesh2014::number_of_te_particles()
{
	return (int) m_te_vortex_particles.size();
}

int mFlow::Ramesh2014::number_of_le_particles()
{
	return (int)m_le_vortex_particles.size();
}

int mFlow::Ramesh2014::number_of_particles()
{
	return number_of_le_particles() + number_of_te_particles();
}


std::vector<double> mFlow::Ramesh2014::get_fourier_terms()
{
	return m_fourier_terms;
}

void mFlow::Ramesh2014::calculate_velocities()
{	// checked.
	// Initialise at zero for sum.
	for (auto & particle_vel : m_te_vortex_particle_velocities) {
		particle_vel.as_array() = { 0, 0 };
	}
	for (auto & particle_vel : m_le_vortex_particle_velocities) {
		particle_vel.as_array() = { 0, 0 };
	}

	// Velocities induced on particles by other particles (N-Body problem)
	if(wake_self_convection){	// checked.
		double vc_size = vortex_core_size();
		for (int i = 0; i < number_of_te_particles(); i++) {
			HBTK::CartesianPoint2D ith_particle_pos = m_te_vortex_particles[i].position;
			for (int j = 0; j < i; j++) {
				HBTK::CartesianVector2D vel;	// Cartesian velocity
				assert(m_te_vortex_particles[j].position != m_te_vortex_particles[i].position);
				vel = unity_vortex_blob_induced_vel(m_te_vortex_particles[j].position,
					ith_particle_pos, vc_size);
				m_te_vortex_particle_velocities[j] += vel * m_te_vortex_particles[i].vorticity;
				m_te_vortex_particle_velocities[i] -= vel * m_te_vortex_particles[j].vorticity;
			}
			for (int j = 0; j < number_of_le_particles(); j++) {
				HBTK::CartesianVector2D vel;	// Cartesian velocity
				vel = unity_vortex_blob_induced_vel(m_le_vortex_particles[j].position,
					ith_particle_pos, vc_size);
				m_le_vortex_particle_velocities[j] += vel * m_te_vortex_particles[i].vorticity;
				m_te_vortex_particle_velocities[i] -= vel * m_le_vortex_particles[j].vorticity;
			}
		}
		for (int i = 0; i < number_of_le_particles(); i++) {
			HBTK::CartesianPoint2D ith_particle_pos = m_le_vortex_particles[i].position;
			for (int j = 0; j < i; j++) {
				HBTK::CartesianVector2D vel;	// Cartesian velocity
				vel = unity_vortex_blob_induced_vel(m_le_vortex_particles[j].position,
					ith_particle_pos, vc_size);
				m_le_vortex_particle_velocities[j] += vel * m_le_vortex_particles[i].vorticity;
				m_le_vortex_particle_velocities[i] -= vel * m_le_vortex_particles[j].vorticity;
			}
		}
	}

	// Velocity induced on particles by the aerofoil's vorticity distribution.
	{
		// We'll integrate over the aerofoil's vorticity distribution using a remapped Gauss quadrature.
		HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
		quad.telles_cubic_remap(-1.0);	// Adapt quadrature for the leading edge singularity.
		double vc_size = vortex_core_size();
		auto u_integrand_outer = [&](double local_foil_pos, HBTK::CartesianPoint2D mes_pnt)->HBTK::CartesianVector2D {
			HBTK::CartesianPoint2D foil_coord = foil_coordinate(local_foil_pos);
			HBTK::CartesianVector2D vel = unity_vortex_blob_induced_vel(mes_pnt, foil_coord, vc_size)
				* vorticity_density(local_foil_pos) * semichord;
			return vel;
		};
#pragma omp parallel for
		for (int i = 0; i < number_of_te_particles(); i++) {
			HBTK::CartesianPoint2D mes_pnt = m_te_vortex_particles[i].position;
			auto integrand_inner = [&](double local_pos) { return u_integrand_outer(local_pos, mes_pnt); };
			m_te_vortex_particle_velocities[i] += quad.integrate(integrand_inner);
		}
#pragma omp parallel for
		for (int i = 0; i < number_of_le_particles(); i++) {
			HBTK::CartesianPoint2D mes_pnt = m_le_vortex_particles[i].position;
			auto integrand_inner = [&](double local_pos) { return u_integrand_outer(local_pos, mes_pnt); };
			m_le_vortex_particle_velocities[i] += quad.integrate(integrand_inner);
		}
	}

	// Velocity induced by free stream.
	for(int i = 0; i < number_of_te_particles(); i++) {
		m_te_vortex_particle_velocities[i] += free_stream_velocity
			+ external_purturbation(m_te_vortex_particles[i].position, time);
	}
	for (int i = 0; i < number_of_le_particles(); i++) {
		m_le_vortex_particle_velocities[i] += free_stream_velocity
			+ external_purturbation(m_le_vortex_particles[i].position, time);
	}

	// We're done. Check that we've done something sane:
	for (auto particle_vel : m_te_vortex_particle_velocities) {
		assert(HBTK::check_finite(particle_vel));
	}
	for (auto particle_vel : m_le_vortex_particle_velocities) {
		assert(HBTK::check_finite(particle_vel));
	}
}

void mFlow::Ramesh2014::convect_particles()
{
	// Forward Euler scheme - justified by Ramesh et al.
	for (int i = 0; i < number_of_te_particles(); i++) {
		m_te_vortex_particles[i].position += m_te_vortex_particle_velocities[i] * delta_t;
	}
	for (int i = 0; i < number_of_le_particles(); i++) {
		m_le_vortex_particles[i].position += m_le_vortex_particle_velocities[i] * delta_t;
	}
	return;
}

void mFlow::Ramesh2014::shed_new_trailing_edge_particle_with_zero_vorticity()
{	// checked.
	VortexGroup2D::Vortex2D particle;
	HBTK::CartesianVector2D velocity;
	if (number_of_te_particles() > 0) {
		particle = m_te_vortex_particles.most_recently_added();
		HBTK::CartesianPoint2D te_coord = foil_coordinate(1);
		// Adjust so that the new particle is 1/3 of the distance
		// between the TE and the last shed particle.
		particle.position -= (2. / 3.) * (particle.position - te_coord);
		particle.vorticity = 0;
	}
	else {
		// Our first particle.
		// Put the particle in the path of foil.
		particle.vorticity = 0;
		particle.position = foil_coordinate(1);
		velocity = -foil_velocity(1);
		velocity += free_stream_velocity + external_purturbation(particle.position, time);
		// 0.5 corresponds to the 1/3 rule used earlier (0.5 + 1 = 3 * 0.5)
		particle.position += velocity * delta_t * 0.5;
	}
	if (!HBTK::check_finite(particle.position)) {
		throw std::runtime_error("Tried to assign infinite position to "
			" new vortex particle. " __FILE__ ": " + std::to_string(__LINE__));
	}
	if (!HBTK::check_finite(particle.position)) {
		throw std::runtime_error("Tried to assign infinite velocity to "
			" new vortex particle. " __FILE__ ": " + std::to_string(__LINE__));
	}
	m_te_vortex_particles.add(particle);
	m_te_vortex_particle_velocities.emplace_back(velocity);
	return;
}

double mFlow::Ramesh2014::total_shed_vorticity()
{	
	return m_wake_vorticity;
}

double mFlow::Ramesh2014::set_total_shed_vorticity_to_wake_vorticity()
{// checked.
	double vorticity = 0;
	for (int i = 0; i < m_te_vortex_particles.size(); i++) {
		vorticity += m_te_vortex_particles[i].vorticity;
	}
	for (int i = 0; i < m_le_vortex_particles.size(); i++) {
		vorticity += m_le_vortex_particles[i].vorticity;
	}
	if (!HBTK::check_finite(vorticity)) {
		std::overflow_error("Infinite vorticity: "
			"The vorticity of the wake was found to be infinite whilst "
			"running Ramesh2014. " __FILE__ ":" + std::to_string(__LINE__));
	}
	m_wake_vorticity = vorticity;
	return vorticity;
}

HBTK::CartesianVector2D mFlow::Ramesh2014::get_particle_induced_velocity(
	HBTK::CartesianPoint2D mes_location)
{	// checked.
	// We want to use reductions which for OpenMP 2.0 means we can't uses classes. Sad.
	double ind_vel_x = 0;
	double ind_vel_y = 0;
	double vc_size = vortex_core_size();
#pragma omp parallel for reduction(+ : ind_vel_x), reduction(+ : ind_vel_y)
	for (int i = 0; i < m_te_vortex_particles.size(); i++) {
		HBTK::CartesianVector2D particle_induced_vel = unity_vortex_blob_induced_vel(
			mes_location,
			m_te_vortex_particles[i].position, vc_size);
		HBTK::CartesianVector2D iv = particle_induced_vel * m_te_vortex_particles[i].vorticity;
		ind_vel_x += iv.x();
		ind_vel_y += iv.y();
	}
#pragma omp parallel for reduction(+ : ind_vel_x), reduction(+ : ind_vel_y)
	for (int i = 0; i < m_le_vortex_particles.size(); i++) {
		HBTK::CartesianVector2D particle_induced_vel = unity_vortex_blob_induced_vel(
			mes_location,
			m_le_vortex_particles[i].position, vc_size);
		HBTK::CartesianVector2D iv = particle_induced_vel * m_le_vortex_particles[i].vorticity;
		ind_vel_x += iv.x();
		ind_vel_y += iv.y();
	}
	HBTK::CartesianVector2D induced_vel({ ind_vel_x, ind_vel_y });
	if (!HBTK::check_finite(induced_vel)) {
		throw std::runtime_error("Trying to evaluate infinite velocity in "
			" flow field. " __FILE__ ": " + std::to_string(__LINE__));
	}
	return induced_vel;
}

double mFlow::Ramesh2014::induced_velocity_normal_to_foil_surface(double local_coordinate)
{	// checked.
	double wash;
	assert(local_coordinate <= 1);
	assert(local_coordinate >= -1);
	if (m_downwash_cache.count(local_coordinate) == 0) {
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(local_coordinate);
		HBTK::CartesianVector2D p_ind_vel = get_particle_induced_velocity(foil_coord);
		HBTK::CartesianVector2D extern_vel = free_stream_velocity + external_purturbation(foil_coord, time);
		double alpha = foil_AoA(time);
		extern_vel.rotate(alpha);
		double alpha_dot = foil_dAoAdt(time);
		double h_dot = foil_dZdt(time);
		double local_camber_slope = camber_slope(local_coordinate);
		p_ind_vel.rotate(alpha);						// Good
		wash = local_camber_slope *
			(extern_vel.x()
				+ h_dot * sin(alpha)					// Good
				+ p_ind_vel.x())						// Good
			- extern_vel.y()
			- alpha_dot * semichord * (local_coordinate - pitch_location)	// Good
			+ h_dot * cos(alpha)						// Good
			- p_ind_vel.y();							// Good
		m_downwash_cache[local_coordinate] = wash;
		if (!HBTK::check_finite(wash)) {
			throw std::runtime_error("Trying to evaluate infinite downwash in "
				" onto aerofoil. " __FILE__ ": " + std::to_string(__LINE__));
		}
	}
	else {
		wash = m_downwash_cache[local_coordinate];
	}
	return wash;
}

void mFlow::Ramesh2014::add_last_shed_te_vortex_to_downwash_cache()
{
	auto & particle = m_te_vortex_particles.most_recently_added();
	for (auto & x : m_downwash_cache) {
		assert(x.first <= 1);
		assert(x.first >= -1);
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(x.first);
		HBTK::CartesianVector2D p_ind_vel = particle.vorticity *
			unity_vortex_blob_induced_vel(foil_coord, particle.position, vortex_core_size());
		double alpha = foil_AoA(time);
		double local_camber_slope = camber_slope(x.first);
		p_ind_vel.rotate(alpha);
		x.second += local_camber_slope * p_ind_vel.x() - p_ind_vel.y();
		if (!HBTK::check_finite(x.second)) {
			throw std::runtime_error("Corrected downwash such that it is "
				" now infinite. Error! " __FILE__ ": " + std::to_string(__LINE__));
		}
	}
}

void mFlow::Ramesh2014::add_last_shed_le_vortex_to_downwash_cache()
{
	auto & particle = m_le_vortex_particles.most_recently_added();
	for (auto & x : m_downwash_cache) {
		assert(x.first <= 1);
		assert(x.first >= -1);
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(x.first);
		HBTK::CartesianVector2D p_ind_vel = particle.vorticity *
			unity_vortex_blob_induced_vel(foil_coord, particle.position, vortex_core_size());
		double alpha = foil_AoA(time);
		double local_camber_slope = camber_slope(x.first);
		p_ind_vel.rotate(alpha);						// Good
		x.second += local_camber_slope * p_ind_vel.x() - p_ind_vel.y();
		if (!HBTK::check_finite(x.second)) {
			throw std::runtime_error("Corrected downwash such that it is "
				" now infinite. Error! " __FILE__ ": " + std::to_string(__LINE__));
		}
	}
}

void mFlow::Ramesh2014::remove_last_shed_te_vortex_from_downwash_cache()
{
	auto & particle = m_te_vortex_particles.most_recently_added();
	for (auto & x : m_downwash_cache) {
		// first is the key - the location at which the downwash
		// function is evaluated on the foil.
		assert(x.first <= 1);
		assert(x.first >= -1);
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(x.first);
		HBTK::CartesianVector2D p_ind_vel = particle.vorticity *
			unity_vortex_blob_induced_vel(foil_coord, particle.position, vortex_core_size());
		double alpha = foil_AoA(time);
		double local_camber_slope = camber_slope(x.first);
		p_ind_vel.rotate(alpha);						// Good
		x.second -= local_camber_slope * p_ind_vel.x() - p_ind_vel.y();
		if (!HBTK::check_finite(x.second)) {
			throw std::runtime_error("Corrected downwash such that it is "
				" now infinite. Error! " __FILE__ ": " + std::to_string(__LINE__));
		}
	}
}


void mFlow::Ramesh2014::adjust_last_shed_vortex_particle_for_kelvin_condition()
{	// checked - possible problems with integration of I_unknown.
	m_te_vortex_particles.most_recently_added().vorticity = 0;
	// The bound vorticity distribution can be decomposed into a part that is a result of
	// the wake for which we have already set the vorticity.
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	auto known_integrand = [&](double theta) -> double {
		double local_x = - cos(theta);		// [0, pi] -> x in [-1,1] -> [LE, TE]
		double T_1 = induced_velocity_normal_to_foil_surface(local_x);
		return T_1 * (cos(theta) - 1) * 2.0 * semichord;
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	double I_known = quad.integrate(known_integrand);
	if (!HBTK::check_finite(I_known)) {
		throw std::runtime_error("Integral I_known"
			"evaluated as infinite. " __FILE__ ": " + std::to_string(__LINE__));
	}

	// And a part that will be cause by our new vortex particle.
	double I_unknown;
	{
		HBTK::CartesianPoint2D particle_pos = m_te_vortex_particles.most_recently_added().position;
		double v_core_size = vortex_core_size();
		auto particle_integrand = [&](double theta) -> double {
			double foil_pos = - cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
			HBTK::CartesianPoint2D foil_coord = foil_coordinate(foil_pos);
			HBTK::CartesianVector2D velocity	// Velocity in foil frame.
				= unity_vortex_blob_induced_vel(foil_coord, particle_pos, v_core_size).rotate(alpha);
			double normal_velocity = camber_slope(foil_pos) * velocity.x() - velocity.y();
			return normal_velocity * (cos(theta) - 1) * 2.0 * semichord;
		};
		// We expect particle_integrand to be singular looking towards the TE.
		quad.telles_quadratic_remap(HBTK::Constants::pi());	
		I_unknown = quad.integrate(particle_integrand);
	}
	if (!HBTK::check_finite(I_unknown)) {
		throw std::runtime_error("Integral I_unknown"
			"evaluated as infinite. " __FILE__ ": " + std::to_string(__LINE__));
	}
	m_te_vortex_particles.most_recently_added().vorticity = - (I_known + total_shed_vorticity()) / (1 + I_unknown);
	m_wake_vorticity += m_te_vortex_particles.most_recently_added().vorticity;
	add_last_shed_te_vortex_to_downwash_cache();
	return;
}

void mFlow::Ramesh2014::shed_new_leading_edge_particle_if_required_and_adjust_vorticities()
{	// CORRECTNESS in free_stream_vel?
	if (!lev_shedding && !shed_zero_strength_levs) { return; }
	// Compute A_0 term.
	auto integrand = [&](double theta) {
		double foil_pos = -cos(theta);
		return induced_velocity_normal_to_foil_surface(foil_pos)
			/ free_stream_velocity.magnitude();
	};
	double A_0 = -0.5 * HBTK::adaptive_simpsons_integrate(integrand, 1.e-7, 0.0, HBTK::Constants::pi()) * 2 / HBTK::Constants::pi();

	if (std::abs(A_0) > critical_leading_edge_suction) {
		// Note that the last TE particle isn't actually removed from simulation.
		double lesp_sign = A_0 > 0 ? 1. : -1.;
		m_wake_vorticity -= m_te_vortex_particles.most_recently_added().vorticity;
		remove_last_shed_te_vortex_from_downwash_cache();
		shed_new_leading_edge_particle_with_zero_vorticity();
		double I_1, I_2, I_3, J_1, J_2, J_3;
		I_1 = known_wake_kelvin_condition_effect_term();
		I_2 = last_tev_kelvin_condition_effect_term();
		I_3 = new_lev_kelvin_condition_effect_term();
		J_1 = known_wake_A_0_effect_term();
		J_2 = last_tev_A_0_effect_term();
		J_3 = new_lev_A_0_effect_term();
		double det = J_3 * (I_2 + 1.) - J_2 * (I_3 + 1.);
		m_te_vortex_particles.most_recently_added().vorticity =
			(1. / det) * (-J_3 * (I_1 + total_shed_vorticity())
				+ (I_3 + 1) * (J_1 - critical_leading_edge_suction * lesp_sign));
		m_le_vortex_particles.most_recently_added().vorticity =
			(1. / det) * (J_2 * (I_1 + total_shed_vorticity())
				- (I_2 + 1) * (J_1 - critical_leading_edge_suction * lesp_sign));
		add_last_shed_le_vortex_to_downwash_cache();
		add_last_shed_te_vortex_to_downwash_cache();
		m_wake_vorticity += m_te_vortex_particles.most_recently_added().vorticity;
		m_wake_vorticity += m_le_vortex_particles.most_recently_added().vorticity;
		if (!HBTK::check_finite(m_wake_vorticity)) {
			throw std::runtime_error("Wake vorticity has blown up! - "
				"evaluated as infinite. " __FILE__ ": " + std::to_string(__LINE__));
		}
		m_shed_lev_at_last_step = true;
	}
	else if (shed_zero_strength_levs) {
		shed_new_leading_edge_particle_with_zero_vorticity();
		m_shed_lev_at_last_step = true;
	}
	else{ m_shed_lev_at_last_step = false; }
	return;
}

void mFlow::Ramesh2014::shed_new_leading_edge_particle_with_zero_vorticity()
{
	VortexGroup2D::Vortex2D particle;
	HBTK::CartesianVector2D velocity;
	if (!m_shed_lev_at_last_step) {
		particle.position = foil_coordinate(-1);
		velocity = -foil_velocity(-1);
		velocity += free_stream_velocity + external_purturbation(particle.position, time);
		// 0.5 corresponds to the 1/3 rule used earlier (0.5 + 1 = 3 * 0.5)
		particle.position += velocity * delta_t * 0.5;
	}
	else {
		particle = m_le_vortex_particles.most_recently_added();
		HBTK::CartesianPoint2D le_coord = foil_coordinate(-1);
		particle.position -= (2. / 3.) * (particle.position - le_coord);
	}
	particle.vorticity = 0;
	m_le_vortex_particles.add(particle);
	m_le_vortex_particle_velocities.emplace_back(velocity);
	return;
}

void mFlow::Ramesh2014::compute_fourier_terms()
{	// Eq.2.3/2.4 Free stream velocity bit correctness?.
	assert(number_of_fourier_terms > 2);
	compute_fourier_terms(0, number_of_fourier_terms);
	return;
}

void mFlow::Ramesh2014::compute_fourier_terms(int minterm, int maxterm) {
	assert(minterm >= 0);
	assert(maxterm > minterm);
	assert(maxterm < number_of_fourier_terms);
	// HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	// quad.linear_remap(0, HBTK::Constants::pi());
	for (int i = minterm; i < maxterm; i++) {
		auto integrand = [&](double theta) {
			double foil_pos = -cos(theta);
			return cos(i * theta) * induced_velocity_normal_to_foil_surface(foil_pos)
				/ free_stream_velocity.magnitude();
		};
		m_fourier_terms[i] = HBTK::adaptive_simpsons_integrate(integrand, 1e-7, 0.0, HBTK::Constants::pi()) * 2. / HBTK::Constants::pi();
	}
	if ((minterm == 0) && (maxterm > 0)) { m_fourier_terms[0] *= -0.5; }
	if (!HBTK::check_finite(m_fourier_terms)) {
		throw std::runtime_error("Fourier terms were evaluated "
			" as infinite! " __FILE__ ": " + std::to_string(__LINE__));
	}
	return;
}

double mFlow::Ramesh2014::bound_vorticity()
{	// Eq.2.7 - free stream velocity correctness?
	return free_stream_velocity.magnitude() *
		semichord * HBTK::Constants::pi() * 
		(2 * m_fourier_terms[0] + m_fourier_terms[1]);
}

double mFlow::Ramesh2014::vorticity_density(double local_pos)
{	// Eq.2.1 - free stream velocity correctness?
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

HBTK::CartesianPoint2D mFlow::Ramesh2014::foil_coordinate(double eta)
{	// checked.
	HBTK::CartesianPoint2D pitch_loc({ pitch_location * semichord, 0 });
	HBTK::CartesianPoint2D static_foil_coord({ eta * semichord, camber_line(eta) });
	HBTK::CartesianPoint2D foil_coord = pitch_loc 
		+ (static_foil_coord - pitch_loc).rotated(-foil_AoA(time));
	foil_coord.y() += foil_Z(time);
	return foil_coord;
}

HBTK::CartesianVector2D mFlow::Ramesh2014::foil_velocity(double eta)
{	// checked.
	double angular_vel = foil_dAoAdt(time);
	HBTK::CartesianVector2D radius = foil_coordinate(eta) - pivot_coordinate();
	HBTK::CartesianVector2D vel;
	vel.y() = - angular_vel * radius.x();
	vel.x() = angular_vel * radius.y();
	vel.y() += foil_dZdt(time);
	return vel;
}

double mFlow::Ramesh2014::aerofoil_leading_edge_suction_force(double density)
{	// free stream vel correctness?
	assert(HBTK::check_finite(m_fourier_terms));
	return pow(free_stream_velocity.magnitude(), 2) * HBTK::Constants::pi() *
		density * 2 * semichord * pow(m_fourier_terms[0], 2);
}

double mFlow::Ramesh2014::aerofoil_normal_force(double density)
{	// free stream vel correctness?
	assert(HBTK::check_finite(m_fourier_terms));
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	std::vector<double>& fourier_derivatives = m_fourier_derivatives;

	double term_1, term_2,
		term_11, term_12, term_121, term_122,
		term_1211, term_1212;

	term_11 = density * HBTK::Constants::pi() * semichord * 2.0 * free_stream_velocity.magnitude();
	term_1211 = free_stream_velocity.rotated(alpha).x() + h_dot * sin(alpha);
	term_1212 = m_fourier_terms[0] + m_fourier_terms[1] / 2;
	term_121  = term_1211 * term_1212;
	term_122 = 2 * semichord * (
		(3. / 4) * fourier_derivatives[0]
		+ (1. / 4) * fourier_derivatives[1]
		+ (1. / 8) * fourier_derivatives[2]);
	term_12 = term_121 + term_122;
	term_1 = term_11 * term_12;
	assert(HBTK::check_finite(term_1));

	// term_2 includes a weakly singular integral. We use the singularity subtraction method.
	double ssm_static = get_particle_induced_velocity(foil_coordinate(-1)).rotate(alpha).x();
	auto integrand = [&](double local_pos)->double {
		double vort_den = vorticity_density(local_pos);
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(local_pos);
		HBTK::CartesianVector2D vel = get_particle_induced_velocity(foil_coord);
		return vort_den * (vel.rotate(alpha).x() - ssm_static);
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	term_2 = quad.integrate(integrand) * semichord * density;
	term_2 += density * ssm_static * free_stream_velocity.magnitude() * 2 * semichord *
		HBTK::Constants::pi() * (m_fourier_terms[0] + m_fourier_terms[1] / 2);
	assert(HBTK::check_finite(term_2));
	return term_1 +term_2;
}

double mFlow::Ramesh2014::aerofoil_moment_about_pitch_location(double density)
{	// NOT CHECKED - SINGULAR INTEGRAL!!!
	assert(HBTK::check_finite(density));
	assert(HBTK::check_finite(m_fourier_terms));
	double alpha = foil_AoA(time);
	double alpha_dot = foil_dAoAdt(time);
	double h_dot = foil_dZdt(time);
	std::vector<double>& fourier_derivatives = m_fourier_derivatives;

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
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(local_pos);
		HBTK::CartesianVector2D vel
			= get_particle_induced_velocity(foil_coord);
		double tangential_velocity = vel.x() * cos(alpha) - vel.y() * sin(alpha);
		return tangential_velocity * vort_den * semichord * (1 + local_pos);
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(40);
	term_3 = quad.integrate(integrand) * semichord * density;
	assert(HBTK::check_finite(term_3));

	return term_1 - term_2 - term_3;
}

std::pair<double, double> mFlow::Ramesh2014::aerofoil_lift_and_drag_coefficients()
{	// free stream vel correctness?
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

void mFlow::Ramesh2014::calc_rate_of_change_of_fourier_terms()
{	// checked.
	assert(delta_t != 0);
	m_fourier_derivatives.resize(m_fourier_terms.size());
	std::vector<double> save_fourier = m_fourier_terms;
	compute_fourier_terms();
	for (int i = 0; i < (int) m_fourier_terms.size(); i++) {
		m_fourier_derivatives[i] = (m_fourier_terms[i] - m_previous_fourier_terms[i]) / delta_t;
	}
	m_fourier_terms = save_fourier;
	assert(HBTK::check_finite(m_fourier_derivatives));
	return;
}

double mFlow::Ramesh2014::known_wake_kelvin_condition_effect_term()
{
	m_le_vortex_particles.most_recently_added().vorticity = 0;
	m_te_vortex_particles.most_recently_added().vorticity = 0;
	// The bound vorticity distribution can be decomposed into a part that is a result of
	// the wake for which we have already set the vorticity.
	auto known_integrand = [&](double theta) -> double {
		double local_x = -cos(theta);		// [0, pi] -> x in [-1,1] -> [LE, TE]
		double T_1 = induced_velocity_normal_to_foil_surface(local_x);
		return T_1 * (cos(theta) - 1) * 2.0 * semichord;
	};	
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	double I_known = quad.integrate(known_integrand); 
	//double I_known = HBTK::adaptive_simpsons_integrate(
	//	known_integrand, 1e-10, 0.0, HBTK::Constants::pi());
	assert(HBTK::check_finite(I_known));
	return I_known;
}

double mFlow::Ramesh2014::last_tev_kelvin_condition_effect_term()
{
	double I_unknown;
	double alpha = foil_AoA(time);
	HBTK::CartesianPoint2D particle_pos = m_te_vortex_particles.most_recently_added().position;
	double v_core_size = vortex_core_size();
	auto particle_integrand = [&](double theta) -> double {
		double foil_pos = -cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(foil_pos);
		HBTK::CartesianVector2D velocity	// Velocity in foil frame.
			= unity_vortex_blob_induced_vel(foil_coord, particle_pos, v_core_size).rotate(alpha);
		double normal_velocity = camber_slope(foil_pos) * velocity.x() - velocity.y();
		return normal_velocity * (cos(theta) - 1) * 2.0 * semichord;
	};
	// We expect particle_integrand to be singular looking towards the TE.
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	quad.telles_quadratic_remap(HBTK::Constants::pi());
	I_unknown = quad.integrate(particle_integrand);
	//I_unknown = HBTK::adaptive_simpsons_integrate(
	//	particle_integrand, 1e-10, 0.0, HBTK::Constants::pi());

	assert(HBTK::check_finite(I_unknown));
	return I_unknown;
}

double mFlow::Ramesh2014::new_lev_kelvin_condition_effect_term()
{
	double I_unknown;
	double alpha = foil_AoA(time);
	HBTK::CartesianPoint2D particle_pos = m_le_vortex_particles.most_recently_added().position;
	double v_core_size = vortex_core_size();
	auto particle_integrand = [&](double theta) -> double {
		double foil_pos = -cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(foil_pos);
		HBTK::CartesianVector2D velocity	// Velocity in foil frame.
			= unity_vortex_blob_induced_vel(foil_coord, particle_pos, v_core_size).rotate(alpha);
		double normal_velocity = camber_slope(foil_pos) * velocity.x() - velocity.y();
		return normal_velocity * (cos(theta) - 1) * 2.0 * semichord;
	};
	// We expect particle_integrand to be singular looking towards the TE.
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	quad.telles_quadratic_remap(0.0);
	I_unknown = quad.integrate(particle_integrand);
	//I_unknown = HBTK::adaptive_simpsons_integrate(
	//	particle_integrand, 1e-10, 0.0, HBTK::Constants::pi());

	assert(HBTK::check_finite(I_unknown));
	return I_unknown;
}

double mFlow::Ramesh2014::known_wake_A_0_effect_term()
{
	m_le_vortex_particles.most_recently_added().vorticity = 0;
	m_te_vortex_particles.most_recently_added().vorticity = 0;
	// The bound vorticity distribution can be decomposed into a part that is a result of
	// the wake for which we have already set the vorticity.
	auto known_integrand = [&](double theta) -> double {
		double local_x = -cos(theta);		// [0, pi] -> x in [-1,1] -> [LE, TE]
		double T_1 = induced_velocity_normal_to_foil_surface(local_x);
		return T_1;
	};
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	double I_known = quad.integrate(known_integrand);
	//double I_known = HBTK::adaptive_simpsons_integrate(
	//	known_integrand, 1e-10, 0.0, HBTK::Constants::pi());
	I_known *= -1. / (HBTK::Constants::pi() * free_stream_velocity.magnitude());
	assert(HBTK::check_finite(I_known));
	return I_known;
}

double mFlow::Ramesh2014::last_tev_A_0_effect_term()
{
	double I_unknown;
	double alpha = foil_AoA(time);
	HBTK::CartesianPoint2D particle_pos = m_te_vortex_particles.most_recently_added().position;
	double v_core_size = vortex_core_size();
	auto particle_integrand = [&](double theta) -> double {
		double foil_pos = -cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(foil_pos);
		HBTK::CartesianVector2D velocity	// Velocity in foil frame.
			= unity_vortex_blob_induced_vel(foil_coord, particle_pos, v_core_size).rotate(alpha);
		double normal_velocity = camber_slope(foil_pos) * velocity.x() - velocity.y();
		return normal_velocity;
	};
	// We expect particle_integrand to be singular looking towards the TE.
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	quad.telles_quadratic_remap(HBTK::Constants::pi());
	I_unknown = quad.integrate(particle_integrand);
	//I_unknown = HBTK::adaptive_simpsons_integrate(
	//	particle_integrand, 1e-10, 0.0, HBTK::Constants::pi());
	I_unknown *= -1. / (HBTK::Constants::pi() * free_stream_velocity.magnitude());

	assert(HBTK::check_finite(I_unknown));
	return I_unknown;
}

double mFlow::Ramesh2014::new_lev_A_0_effect_term()
{
	double I_unknown;
	double alpha = foil_AoA(time);
	HBTK::CartesianPoint2D particle_pos = m_le_vortex_particles.most_recently_added().position;
	double v_core_size = vortex_core_size();
	auto particle_integrand = [&](double theta) -> double {
		double foil_pos = -cos(theta);	// Foil local coordinate in -1-> LE, 1->TE
		HBTK::CartesianPoint2D foil_coord = foil_coordinate(foil_pos);
		HBTK::CartesianVector2D velocity	// Velocity in foil frame.
			= unity_vortex_blob_induced_vel(foil_coord, particle_pos, v_core_size).rotate(alpha);
		double normal_velocity = camber_slope(foil_pos) * velocity.x() - velocity.y();
		return normal_velocity;
	};
	// We expect particle_integrand to be singular looking towards the TE.
	HBTK::StaticQuadrature quad = HBTK::gauss_legendre(50);
	quad.linear_remap(0, HBTK::Constants::pi());
	quad.telles_quadratic_remap(0.0);
	I_unknown = quad.integrate(particle_integrand);
	//I_unknown = HBTK::adaptive_simpsons_integrate(
	//	particle_integrand, 1e-10, 0.0, HBTK::Constants::pi());
	I_unknown *= -1. / (HBTK::Constants::pi() * free_stream_velocity.magnitude());

	assert(HBTK::check_finite(I_unknown));
	return I_unknown;
}

double mFlow::Ramesh2014::vortex_core_size() const
{	// checked.
	double delta_t_star_times_c = delta_t * free_stream_velocity.magnitude();
	return delta_t_star_times_c * 1.3;
}

HBTK::CartesianPoint2D mFlow::Ramesh2014::pivot_coordinate()
{	// checked.
	assert(std::abs(pitch_location) <= 1);
	HBTK::CartesianPoint2D loc;
	loc.x() = pitch_location * semichord;
	loc.y() = foil_Z(time);
	return loc;
}

HBTK::CartesianVector2D mFlow::Ramesh2014::unity_vortex_blob_induced_vel(
	const HBTK::CartesianPoint2D & mes_point,
	const HBTK::CartesianPoint2D & vortex_point,
	double vortex_size)
{	// checked.
	assert(mes_point != vortex_point);
	double denominator;
	HBTK::CartesianVector2D diff = mes_point - vortex_point;
	denominator = sqrt(pow(pow(diff.x(), 2) + pow(diff.y(), 2), 2)
		+ pow(vortex_size, 4)) * 2. * HBTK::Constants::pi();
	HBTK::CartesianVector2D vel;
	vel.x() = diff.y() / denominator;
	vel.y() = -diff.x() / denominator;
	return vel;
}
