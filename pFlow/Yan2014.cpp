#include "stdafx.h"
#include "Yan2014.h"

#include "HBTK/ConformalMapping.h"
#include "HBTK/Constants.h"
#include "HBTK/PotentialFlowDistributions.h"

mFlow::Yan2014::Yan2014()
{
}


mFlow::Yan2014::~Yan2014()
{
}

double mFlow::Yan2014::particle_arg_eta_plain(const vortex_particle & particle)
{
	return std::arg(particle.x);
}

void mFlow::Yan2014::compute_eta_plane_velocities(double time)
{
	int n_particles = (int) m_vortex_particles.size();
	std::vector<double> rad, arg, vrad, vtan;
	rad.resize(n_particles);
	arg.resize(n_particles);
	vrad.resize(n_particles);
	vtan.resize(n_particles);
	for (int i = 0; i < n_particles; i++) {
		rad[i] = abs(m_vortex_particles[i].x);
		arg[i] = std::arg(m_vortex_particles[i].x);
		vtan[i] = 0;
		vrad[i] = 0;
	}

	// p* is arg in eta plane from centre of disk
	// r* is radius in eta plane from centre of disk
	// G* is vorticity
	// j,k are particle we're examining the induced velocity on
	// and inducing particle repectively.
	auto kernel = [=](double pj, double pk, double rj, double rk, double Gk)->std::pair<double, double> {
		double qr, qt; // radial and tangential velocities.
		double b2 = pow(semichord, 2); // Semispan squared.
		double den1 = rj * rj - 2 * rj * b2 * cos(pj - pk) / (4 * rj) + pow(b2 / (4 * rk), 2);
		double den2 = rj * rj - 2 * rj * rk * cos(pj - pk) + rk * rk;
		// Eq30a
		double rterm_1 = -(Gk * sin(pj - pk)) / (2 * HBTK::Constants::pi());
		double rterm_2 = rj / den2 - b2 / (4 * rk * den1);
		// Eq30b
		double tterm_1 = -Gk / (4 * HBTK::Constants::pi() * rj);
		double tterm_21 = (rj*rj - pow(b2 / (4 * rk), 2)) / den1;
		double tterm_22 = (rj * rj - rk * rk) / den1;

		return std::make_pair(rterm_1 * rterm_2, tterm_1 * (tterm_21 - tterm_22));
	};

	for (int i = 0; i < n_particles; i++) {
		for (int j = 0; j < i; j++) {
			double pi, pj, ri, rj, Gi, Gj, vr, vt;
			pi = std::arg(m_vortex_particles[i].x);
			ri = std::arg(m_vortex_particles[j].x);
			rj = abs(m_vortex_particles[i].x);
			rj = abs(m_vortex_particles[j].x);
			Gi = m_vortex_particles[i].vorticity;
			Gj = m_vortex_particles[j].vorticity;
			std::tie(vr, vt) = kernel(pi, pj, ri, rj, Gj);
			vrad[i] += vr;
			vtan[i] += vt;
			std::tie(vr, vt) = kernel(pj, pi, rj, ri, Gi);
			vrad[j] += vr;
			vtan[j] += vt;
		}
	}

	// Now the non-circulatory components.
	for (int i = 0; i < n_particles; i++) {
		double pj, rj, b2, alpha, alpha_dot, h_dot, U;
		U = free_stream_velocity;
		pj = std::arg(m_vortex_particles[i].x);
		rj = abs(m_vortex_particles[i].x);
		b2 = pow(semichord, 2);
		alpha = foil_AoA(time);
		alpha_dot = foil_dAoAdt(time);
		h_dot = foil_dZdt(time);

		double t1 = b2 * (U * sin(alpha) + h_dot * cos(alpha) + 
			pitch_location * semichord * b2) / (2 * rj * rj);
		double t2 = b2 * b2 * alpha_dot / (8 * rj * rj * rj);
		vrad[i] += -t1 * sin(pj) - t2 * sin(2 * pj);
		vtan[i] += t1 * cos(pj) + t2 * cos(2 * pj);
	}

	// Now put these velocities into the global coordinate system.
	for (int i = 0; i < n_particles; i++) {
		double pj, b24rj2, rj, qrj, qtj;
		pj = std::arg(m_vortex_particles[i].x);
		rj = abs(m_vortex_particles[i].x);
		b24rj2 = pow(semichord, 2) / (4 * rj * rj);
		qrj = vrad[i];
		qtj = vtan[i];
		double x_dot, z_dot;
		double den = 1. - b24rj2 * cos(pj) + pow(b24rj2, 2);
		x_dot = (qrj * (1 - b24rj2) * cos(pj) - qtj * (1 + b24rj2) * sin(pj)) / den;
		z_dot = (qrj * (1 + b24rj2) * sin(pj) + qtj * (1 - b24rj2) * cos(pj)) / den;

		double X_dot, Z_dot, alpha;
		alpha = foil_AoA(time);
		X_dot = x_dot * cos(alpha) + z_dot * sin(alpha);
		Z_dot = -x_dot * sin(alpha) + z_dot * cos(alpha);

		m_vortex_particles[i].vX = X_dot + free_stream_velocity + 
			pitch_location * semichord * foil_dAoAdt(time) * sin(alpha);
		m_vortex_particles[i].vZ = Z_dot + foil_dZdt(time) +
			pitch_location * semichord * foil_dAoAdt(time) * cos(alpha);
	}

	return;
}

void mFlow::Yan2014::update_eta_plane_vortex_particle_coordinates(double time)
{
	double ca, sa, h;
	ca = cos(foil_AoA(time));
	sa = sin(foil_AoA(time));
	h = foil_Z(time);
	for (auto & particle : m_vortex_particles) {
		double local_X, local_Z;
		local_X = particle.X * ca + (particle.Z - h) * sa;
		local_Z = particle.X * -sa + (particle.Z - h) * ca;
		std::complex<double> X = local_X + HBTK::Constants::i() * local_Z;
		particle.x = X / 2.0 + sqrt(X * X / 4. - 1.); // We can do this because we've a plate.
	}
	return;
}
