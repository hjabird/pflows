#pragma once
/*////////////////////////////////////////////////////////////////////////////
Sclavounos1987.cpp

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

#include <cassert>
#include <complex>
#include <functional>
#include <vector>

#include "HBTK/StructuredMeshBlock3D.h"
#include "HBTK/StructuredValueBlockND.h"

namespace mFlow {
	class Yan2014
	{
	public:
		Yan2014();
		~Yan2014();

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

		// A representation of the vortex particles
		struct vortex_particle {
			double X, Z;			// Global cartesian coordinate
			std::complex<double> x; // Foil eta plane coordinate
			double vX, vZ;			// Global velocities.
			double vorticity;		// vorticity of the particle.
		};
		std::vector<vortex_particle> m_vortex_particles;	// Vortex particle array. 
		int number_of_particles() const;					// Number of particles

		HBTK::StructuredMeshBlock3D get_foil_location();
		HBTK::StructuredMeshBlock3D get_3d_vortex_locations();
		HBTK::StructuredValueBlockND<3, double> get_vorticities();

	private:
		// Compute the velocities.
		void compute_eta_plane_velocities();

		// Add a new particle.
		void add_new_vortex_particle();

		void update_eta_plane_vortex_particle_coordinates();
		void forward_euler_convection();

		// Coordinate mappings.
		std::complex<double> map_inertial_to_foil(double X, double Y);
		std::complex<double> map_foil_to_eta(std::complex<double> X);
		std::complex<double> map_eta_to_foil(std::complex<double> eta);
		std::complex<double> map_foil_to_real(std::complex<double> eta);
	};
} // End mFlow

