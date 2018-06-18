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
#include <unordered_map>
#include <vector>

#include <HBTK/CartesianPoint.h>
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
		VortexGroup2D m_te_vortex_particles;
		std::vector<HBTK::CartesianVector2D> m_te_vortex_particle_velocities;
		VortexGroup2D m_le_vortex_particles;
		std::vector<HBTK::CartesianVector2D> m_le_vortex_particle_velocities;
		// Returns number of particles in wake.
		int number_of_te_particles();
		int number_of_le_particles();
		int number_of_particles();


		// Number of fourier terms to represent aerofoil vorticity distribution.
		int number_of_fourier_terms;
		
		// For -1 = LE, 1 = TE, get coordinate of point on foil.
		HBTK::CartesianPoint2D foil_coordinate(double eta);

		// For -1 = LE, 1 = TE, get the velocity of a point on the foil described
		// by foil_coordinate(eta)
		HBTK::CartesianVector2D foil_velocity(double eta);

		// Shed leading edge vortices
		bool lev_shedding;
		// The parameter that governs the onset of shedding.
		double critical_leading_edge_suction;

		// Force:

		// Compute the Blasius leading edge suction force. Eq2.31
		double aerofoil_leading_edge_suction_force(double density);
		// Compute the normal force on the aerofoil. Eq2.30
		double aerofoil_normal_force(double density);
		// Compute the moment about the pitch location Eq2.32
		double aerofoil_moment_about_pitch_location(double density);	// POSSIBLY BAD.
		// Compute the lift and drag coefficients.
		// Std::pair<(lift coeff), (drag coeff)> is returned.
		std::pair<double, double> aerofoil_lift_and_drag_coefficients();

		// Detail:

		// The fourier terms curently representing the vorticity distribution
		// over the aerofoil.
		std::vector<double> get_fourier_terms();

		// Get the total vorticity that has been shed into the wake.
		double total_shed_vorticity();
		// Reset vorticity based on the sum of the vorticities in the particle collector.
		double set_total_shed_vorticity_to_wake_vorticity();

		// Compute the aerofoil's bound vorticity.
		double bound_vorticity();

	private:
		// Compute the velocities of the vortex particles in the wake.
		void calculate_velocities();

		// Convect the wake vortex particles using previosly calculated velocities
		// (effectively a forward Euler scheme.)
		void convect_particles();

		// Shed a new vortex particle from the trailing edge. Vorticity not calculated.
		void shed_new_trailing_edge_particle_with_zero_vorticity();

		// Get the velocity induced at a point by the all the vortex particles
		// in the wake. Invalid if location is that of a particle (singular).
		HBTK::CartesianVector2D get_particle_induced_velocity(
			HBTK::CartesianPoint2D mesurement_location);

		// Get induced velocity normal to foil surface - Eq.2.5
		double induced_velocity_normal_to_foil_surface(double local_coordinate);
		std::unordered_map<double, double> m_downwash_cache;
		void add_last_shed_te_vortex_to_downwash_cache();
		void add_last_shed_le_vortex_to_downwash_cache();
		void remove_last_shed_te_vortex_from_downwash_cache();

		// Adjust the vorticity of the last shed particle to enforce 
		// Kelvin's condition Eq2.6
		void adjust_last_shed_vortex_particle_for_kelvin_condition();

		// If we're shedding LEVs and we've hit the criticality condition, 
		// shed a LEV, and fiddle with the TEV and LEV vorticities to 
		// satify A_0 = LESPcrit and kelvin condition.
		void shed_new_leading_edge_particle_if_required_and_adjust_vorticities();
		void shed_new_leading_edge_particle_with_zero_vorticity();

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

		// Compute the terms used for applying the Kelvin condition.
		double known_wake_kelvin_condition_effect_term(); // I_1
		double last_tev_kelvin_condition_effect_term();	// I_2
		double new_lev_kelvin_condition_effect_term(); // I_3
		double known_wake_A_0_effect_term(); // J_1
		double last_tev_A_0_effect_term(); // J_2
		double new_lev_A_0_effect_term(); // J_3

		// Wake vorticity
		double m_wake_vorticity;

		// The size of the vortex core used for by Ramesh. Eq2.14 and Eq2.15
		double vortex_core_size() const;

		// Get the coordinate for the pivot at any point in time.
		HBTK::CartesianPoint2D pivot_coordinate( void );

		// The velocity induced at point (x_mes, y_mes) by a vortex blob at (x_vor, y_vor)
		// with unity total vorticity. Eq2.12 and Eq2.13
		HBTK::CartesianVector2D unity_vortex_blob_induced_vel(
			const HBTK::CartesianPoint2D & mes_point, 
			const HBTK::CartesianPoint2D & vortex_point,
			double vortex_core_size);
	};
} // End namespace mFlow

