#pragma once
/*////////////////////////////////////////////////////////////////////////////
Ramesh2014.h

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

#include <functional>
#include <tuple>
#include <vector>

#include <HBTK/CartesianVector.h>

#include "VortexGroup2D.h"

namespace mFlow {
	class Ramesh2014
	{
	public:

		Ramesh2014();
		~Ramesh2014();

		HBTK::CartesianVector2D free_stream_velocity;	// (no unit)
		double semichord;				// Half the chord length.
		double pitch_location;			// in [-1, 1]. -1->LE, 1->TE
		// Kinematics.	Functions with respect to time.
		std::function<double(double)> foil_Z;		// h(t)
		std::function<double(double)> foil_dZdt;	// h_dot(t)
		std::function<double(double)> foil_AoA;		// alpha(t) (radians)
		std::function<double(double)> foil_dAoAdt;	// alpha_dot(t) (rad / time)

		// Foil geometry:
		std::function<double(double)> camber_line;	// Defined in [-1, 1]
		std::function<double(double)> camber_slope;	// Defined in [-1, 1]

		// Should the wake convect itself? Deafault is yes.
		bool wake_self_convection;

		double time;				// Current simulation time value.
		double delta_t;				// time = time + delta_t for time stepping.
		// Initialises data structs. Throws exception for bad input values:
		// Exceptions: std::domain_error with explanation message.
		void initialise();			
		// Advance one step in the simulation.
		// Sheds particle, calcs vorticity, convects everything, increments time.
		void advance_one_step();	

		// Vortex Particles:
		VortexGroup2D m_vortex_particles;
		std::vector<HBTK::CartesianVector2D> m_vortex_particle_velocities;
		// Returns number of particles in wake.
		int number_of_particles();


		// Number of fourier terms to represent aerofoil vorticity distribution.
		int number_of_fourier_terms;
		
		// For -1 = LE, 1 = TE, get coordinate of point on foil.
		std::pair<double, double> foil_coordinate(double eta);

		// For -1 = LE, 1 = TE, get the velocity of a point on the foil described
		// by foil_coordinate(eta)
		std::pair<double, double> foil_velocity(double eta);

		// Force:

		// Compute the Blasius leading edge suction force. Eq2.31
		double aerofoil_leading_edge_suction_force(double density);
		// Compute the normal force on the aerofoil. Eq2.30
		double aerofoil_normal_force(double density);
		// Compute the moment about the pitch location Eq2.32
		double aerofoil_moment_about_pitch_location(double density);
		// Compute the lift and drag coefficients.
		// Std::pair<(lift coeff), (drag coeff)> is returned.
		std::pair<double, double> aerofoil_lift_and_drag_coefficients();

		// Detail:

		// The fourier terms curently representing the vorticity distribution
		// over the aerofoil.
		std::vector<double> get_fourier_terms();

		// Get the total vorticity of all vortex particles in the wake.
		double total_shed_vorticity();

		// Compute the aerofoil's bound vorticity.
		double bound_vorticity();

	private:
		// Compute the velocities of the vortex particles in the wake.
		void calculate_velocities();

		// Convect the wake vortex particles using previosly calculated velocities
		// (effectively a forward Euler scheme.)
		void convect_particles();

		// Shed a new vortex particle from the trailing edge. Vorticity not calculated.
		void shed_new_particle();


		// Get the velocity induced at a point by the all the vortex particles
		// in the wake. Invalid if location is that of a particle (singular).
		std::pair<double, double> get_particle_induced_velocity(double x, double y);

		// Adjust the vorticity of the last shed particle to enforce 
		// Kelvin's condition Eq2.6
		void adjust_last_shed_vortex_particle_for_kelvin_condition();

		// For a wake, compute the fourier terms describing the aerofoils
		// vorticity distribution.
		void compute_fourier_terms();


		// For local pos in [-1,1] with -1->LE, 1->TE, evaluate the foil's vorticity
		// density function. Do not evaluate at leading edge (singular).
		double vorticity_density(double local_pos);

		// Vector to hold the fourier terms that describe the vorticity 
		// over the aerofoil as described by Eq2.1 and Eq2.2
		std::vector<double> m_fourier_terms;
		// And the terms from the step before incase we want to compute the 
		// time derivative.
		std::vector<double> m_previous_fourier_terms;
		// A method to extract the time rate of change of fourier terms
		std::vector<double> rate_of_change_of_fourier_terms();


		// The size of the vortex core used for by Ramesh. Eq2.14 and Eq2.15
		double vortex_core_size() const;

		// The velocity induced at point (x_mes, y_mes) by a vortex blob at (x_vor, y_vor)
		// with unity total vorticity. Eq2.12 and Eq2.13
		std::pair<double, double> unity_vortex_blob_induced_vel(
			double x_mes, double y_mes, double x_vor, double y_vor, double vortex_core_size);
	};
} // End namespace mFlow

