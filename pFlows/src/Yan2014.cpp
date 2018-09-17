#include "Yan2014.h"
/*////////////////////////////////////////////////////////////////////////////
Yan2014.cpp

Method from the paper "Geometrically-exact unsteady model for airfoils
undergoing large amplitude maneuvers" Aerospace Science and Technology
2014.
https://doi.org/10.1016/j.ast.2014.09.021

Copyright 2017 HJA Bird

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

#include <iostream>

#include "HBTK/Checks.h"
#include "HBTK/ConformalMapping.h"
#include "HBTK/Constants.h"
#include "HBTK/PotentialFlowDistributions.h"

mFlow::Yan2014::Yan2014()
	: time(0),
	delta_t(0)
{
}


mFlow::Yan2014::~Yan2014()
{
}

void mFlow::Yan2014::advance_one_step()
{
	assert(delta_t != 0);
	time += delta_t;
	assert(HBTK::check_finite(time));
	assert(foil_AoA);
	assert(foil_dAoAdt);
	assert(foil_Z);
	assert(foil_dZdt);
	assert(HBTK::check_finite(pitch_location));
	update_eta_plane_vortex_particle_coordinates();
	add_new_vortex_particle();
	compute_eta_plane_velocities();
	forward_euler_convection();
}


int mFlow::Yan2014::number_of_particles() const
{
	return (int) m_vortex_particles.size();
}

HBTK::StructuredMeshBlock3D mFlow::Yan2014::get_foil_location()
{
	HBTK::StructuredMeshBlock3D mesh;
	mesh.set_extent({ 2, 1, 1 });
	double TeX, TeZ, LeX, LeZ;
	LeX = - cos(foil_AoA(time)) * semichord;
	TeX = -LeX;
	TeZ = foil_Z(time) + sin(foil_AoA(time)) * semichord;
	LeZ = foil_Z(time) - sin(foil_AoA(time)) * semichord;

	mesh.set_coord({ 0, 0, 0 }, { TeX, 0.0, TeZ });
	mesh.set_coord({ 1, 0, 0 }, { LeX, 0.0, TeZ });

	return HBTK::StructuredMeshBlock3D();
}

HBTK::StructuredMeshBlock3D mFlow::Yan2014::get_3d_vortex_locations()
{
	HBTK::StructuredMeshBlock3D mesh;
	int n_particles = (int)m_vortex_particles.size();
	mesh.set_extent({ n_particles, 1, 1 });
	for (int i = 0; i < n_particles; i++) {
		vortex_particle & particle = m_vortex_particles[i];
		mesh.set_coord({ i, 0, 0 }, { particle.X, 0.0, particle.Z });
	}
	return mesh;
}

HBTK::StructuredValueBlockND<3, double> mFlow::Yan2014::get_vorticities()
{
	HBTK::StructuredValueBlockND<3, double> mesh;
	int n_particles = (int)m_vortex_particles.size();
	mesh.extent({ n_particles, 1, 1 });
	for (int i = 0; i < n_particles; i++) {
		vortex_particle & particle = m_vortex_particles[i];
		mesh.value({ i, 0, 0 }) = particle.vorticity;
	}
	return mesh;
}

void mFlow::Yan2014::compute_eta_plane_velocities()
{
	int n_particles = (int) m_vortex_particles.size();

	std::cout << "Particles: (" << n_particles << ")\n";
	for (int i = 0; i < n_particles; i++) {
		std::cout << m_vortex_particles[i].x << "\n";
	}
	std::cout << "\n";

	// Initialise this point in time's data.
	std::vector<double> rad, arg, vrad, vtan;
	rad.resize(n_particles);	// Radius of particle
	arg.resize(n_particles);	// arg (angle) of particle
	vrad.resize(n_particles);	// Radial velocity of particle
	vtan.resize(n_particles);	// Tangential velocity of particle
	for (int i = 0; i < n_particles; i++) {
		rad[i] = std::abs(m_vortex_particles[i].x);
		assert(rad[i] > semichord / 2);
		arg[i] = std::arg(m_vortex_particles[i].x);
		vtan[i] = 0;
		vrad[i] = 0;
	}
	assert(HBTK::check_finite(rad));
	assert(HBTK::check_finite(arg));

	// Compute the vorticity needed to satisfy Kutta condition and apply
	// to particle that was last convected from the trailing edge.
	{
		double accumulator = 0.0;
		double total_vorticity_acc = 0;
		for (int i = 0; i < n_particles - 1; i++) {
			double rj = rad[i];
			double pj = arg[i];
			double b = semichord;
			double t1 = m_vortex_particles[i].vorticity / (2 * HBTK::Constants::pi() * b);
			total_vorticity_acc += m_vortex_particles[i].vorticity;
			double t2n = rj * rj - b * b / 4;
			double t2d = rj * rj - rj * b * cos(pj) + b * b / 4;
			double vortex_contribution = t1 * t2n / t2d;
			accumulator += vortex_contribution;
		}
		accumulator += free_stream_velocity * sin(foil_AoA(time))
			+ foil_dZdt(time) * cos(foil_AoA(time))
			+ (0.5 + pitch_location) * semichord * foil_dAoAdt(time);
		assert(m_vortex_particles.back().vorticity == 0);
		double rn = rad.back();
		double pn = arg.back();
		double mult = 2. * HBTK::Constants::pi() * semichord *
			((rn * rn - rn * semichord * cos(pn) + pow(semichord / 2., 2)) /
				rn * rn - pow(semichord / 2, 2));
		m_vortex_particles.back().vorticity = accumulator * mult;
		total_vorticity_acc -= accumulator * mult;
		std::cout << "System vorticity: " << total_vorticity_acc << ".\n";
	}

	// The circulatory influence on the vortex particles.
	{
		// Inter-particle velocity influence (circulatory part of problem)
		// p* is arg in eta plane from centre of disk
		// r* is radius in eta plane from centre of disk
		// G* is vorticity
		// j,k are particle we're examining the induced velocity on
		// and inducing particle repectively.
		auto kernel = [=](double pj, double pk, double rj, double rk, double Gk)->std::pair<double, double> {
			double b2 = pow(semichord, 2); // Semispan squared.
			double den1 = pow(rj * rj - 2 * rj * b2 * cos(pj - pk) / (4 * rk) + pow(b2 / (4 * rk), 2), 2);
			double den2 = pow(rj * rj - 2 * rj * rk * cos(pj - pk) + rk * rk, 2);
			// Eq30a
			double rterm_1 = -(Gk * sin(pj - pk)) / (2 * HBTK::Constants::pi());
			double rterm_2 = rj / den2 - (b2 / (4 * rk)) / den1;
			// Eq30b
			double tterm_1 = -Gk / (4 * HBTK::Constants::pi() * rj);
			double tterm_21 = (rj*rj - pow(b2 / (4 * rk), 2)) / den1;
			double tterm_22 = (rj * rj - rk * rk) / den1;

			return std::make_pair(rterm_1 * rterm_2, tterm_1 * (tterm_21 - tterm_22));
		};

		// Compute the convective effects on particles.
		for (int i = 0; i < n_particles; i++) {
			double pi, ri, Gi;
			pi = arg[i];
			ri = rad[i];
			Gi = m_vortex_particles[i].vorticity;
			for (int j = 0; j < i; j++) {
				double pj, rj, Gj, vr, vt;
				pj = arg[j];
				rj = rad[j];
				Gj = m_vortex_particles[j].vorticity;
				std::tie(vr, vt) = kernel(pi, pj, ri, rj, Gj);
				vrad[i] += vr;
				vtan[i] += vt;
				std::tie(vr, vt) = kernel(pj, pi, rj, ri, Gi);
				vrad[j] += vr;
				vtan[j] += vt;

			}
			// Additionally, the particles bound vorticity on its wake
			// vortex pair:
			// vtan[i] -= 2 * ri * Gi / (HBTK::Constants::pi() * (4 * ri * ri - pow(semichord, 2)));
		}
		assert(HBTK::check_finite(vrad));
		assert(HBTK::check_finite(vtan));
	}

	// Now the non-circulatory components.
	{
		double alpha = foil_AoA(time);
		double alpha_dot = foil_dAoAdt(time);
		double h_dot = foil_dZdt(time);
		double U = free_stream_velocity;

		for (int i = 0; i < n_particles; i++) {
			double pj, rj, b2;
			pj = arg[i];
			rj = rad[i];
			b2 = pow(semichord, 2);
			double t1 = b2 * (U * sin(alpha) + h_dot * cos(alpha) +
				pitch_location * semichord * alpha_dot) / (2 * rj * rj);
			double t2 = b2 * b2 * alpha_dot / (8 * rj * rj * rj);
			vrad[i] += -t1 * sin(pj) - t2 * sin(2 * pj);
			vtan[i] += t1 * cos(pj) + t2 * cos(2 * pj);
		}
	}

	// Now put these velocities into the global coordinate system.
	for (int i = 0; i < n_particles; i++) {
		double pj, b24rj2, rj, qrj, qtj;
		pj = arg[i];
		rj = rad[i];
		b24rj2 = pow(semichord, 2) / (4 * rj * rj);
		qrj = vrad[i];
		qtj = vtan[i];
		double x_dot, z_dot;
		double den = 1. - 2 * b24rj2 * cos(pj) + pow(b24rj2, 2);
		x_dot = (qrj * (1 - b24rj2) * cos(pj) - qtj * (1 + b24rj2) * sin(pj)) / den;
		z_dot = (qrj * (1 + b24rj2) * sin(pj) + qtj * (1 - b24rj2) * cos(pj)) / den;

		double X_dot, Z_dot, alpha;
		alpha = foil_AoA(time);
		X_dot = x_dot * cos(alpha) + z_dot * sin(alpha);
		Z_dot = -x_dot * sin(alpha) + z_dot * cos(alpha);
		assert(HBTK::check_finite(X_dot));
		assert(HBTK::check_finite(Z_dot));

		m_vortex_particles[i].vX = X_dot + free_stream_velocity + 
			pitch_location * semichord * foil_dAoAdt(time) * sin(alpha);
		m_vortex_particles[i].vZ = Z_dot + foil_dZdt(time) +
			pitch_location * semichord * foil_dAoAdt(time) * cos(alpha);
	}

	return;
}

void mFlow::Yan2014::add_new_vortex_particle()
{

	vortex_particle te_particle, last_shed_particle;
	if (m_vortex_particles.size() == 0) {
		std::complex<double> pos_guess = semichord + delta_t * free_stream_velocity * 0.5;
		pos_guess = map_foil_to_real(pos_guess);
		last_shed_particle.X = pos_guess.real();
		last_shed_particle.Z = pos_guess.imag();
	}
	else {
		last_shed_particle = m_vortex_particles.back();
	}


	std::complex<double> te_pos_eta = semichord / 2;
	std::complex<double> pos = map_eta_to_foil(te_pos_eta);
	pos = map_foil_to_real(pos);
	pos += (1. / 3.) * (last_shed_particle.X - pos.real());
	pos += HBTK::Constants::i() *(1. / 3.) * (last_shed_particle.Z - pos.imag());
	te_particle.X = pos.real();
	te_particle.Z = pos.imag();
	pos = map_inertial_to_foil(pos.real(), pos.imag());
	pos = map_foil_to_eta(pos);
	te_particle.x = pos;
	te_particle.vorticity = 0;
	te_particle.vX = 0;
	te_particle.vZ = 0;
	m_vortex_particles.emplace_back(te_particle);
}

void mFlow::Yan2014::update_eta_plane_vortex_particle_coordinates()
{
	double ca, sa, h;
	ca = cos(foil_AoA(time));
	sa = sin(foil_AoA(time));
	h = foil_Z(time);
	for (auto & particle : m_vortex_particles) {
		std::complex<double> x = map_inertial_to_foil(particle.X, particle.Z);
		particle.x = map_foil_to_eta(x);
	}
	return;
}

void mFlow::Yan2014::forward_euler_convection()
{
	assert(HBTK::check_finite(delta_t));
	for (auto &particle : m_vortex_particles) {
		assert(HBTK::check_finite(particle.vX));
		assert(HBTK::check_finite(particle.vZ));
		particle.X += particle.vX * delta_t;
		particle.Z += particle.vZ * delta_t;
	}
}

std::complex<double> mFlow::Yan2014::map_inertial_to_foil(double X, double Y)
{
	std::complex<double> foil_coord = 0;
	foil_coord = HBTK::Constants::i() * (Y - foil_Z(time));
	foil_coord += X;
	foil_coord *= exp(HBTK::Constants::i() * foil_AoA(time));
	return foil_coord;
}

std::complex<double> mFlow::Yan2014::map_foil_to_eta(std::complex<double> X)
{
	std::complex<double> eta = (X + sqrt(X*X - semichord * semichord)) / 2.;
	return eta;
}

std::complex<double> mFlow::Yan2014::map_eta_to_foil(std::complex<double> eta)
{
	return eta + semichord * semichord / (4.0 * eta);
}

std::complex<double> mFlow::Yan2014::map_foil_to_real(std::complex<double> x)
{
	x *= exp(- HBTK::Constants::i() * foil_AoA(time));
	x += HBTK::Constants::i() * foil_Z(time);
	return x;
}
