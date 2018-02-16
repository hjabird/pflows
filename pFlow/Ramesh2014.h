#pragma once

#include <functional>
#include <tuple>
#include <vector>

namespace mFlow {
	class Ramesh2014
	{
	public:
		/*
		Right now I'm using HBTK::PointVortex.
		This is a vortex singularity, not a vortex blob 
		as used in Ramesh et al.
		FYI.
		*/


		Ramesh2014();
		~Ramesh2014();

		// Kinematics.
		double free_stream_velocity;
		double semichord;
		double pitch_location; // in [-1, 1]. -1->LE, 1->TE
		std::function<double(double)> foil_Z;
		std::function<double(double)> foil_dZdt;
		std::function<double(double)> foil_AoA;
		std::function<double(double)> foil_dAoAdt;

		double time;
		double delta_t;
		void advance_one_step();

		// vortex particles
		struct vortex_particle {
			double x, y;
			double vx, vy;
			double vorticity;
		};
		std::vector<vortex_particle> m_vortex_particles;
		int num_particles();

		// Fourier series settings.
		int num_fourier_terms;
		
	private:
		void calculate_velocities();
		void convect_particles();
		void shed_new_particle();
		double shed_vorticity();
		std::pair<double, double> get_particle_induced_velocity(double x, double y);
		void adjust_last_shed_vortex_particle_for_kelvin_condition();

		std::vector<double> m_fourier_terms;
		void compute_fourier_terms();
		double bound_vorticity();

		// For -1 = LE, 1 = TE, get coordinate of point on foil.
		std::pair<double, double> foil_coordinate(double eta);
		std::pair<double, double> foil_velocity(double eta);
	};
} // End namespace mFlow

