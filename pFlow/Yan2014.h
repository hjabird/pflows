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
		std::function<double(double)> foil_X;
		std::function<double(double)> foil_dXdt;
		std::function<double(double)> foil_Z;
		std::function<double(double)> foil_dZdt;
		std::function<double(double)> foil_AoA;
		std::function<double(double)> foil_dAoAdt;
		
		// A representation of the vortex particles
		struct vortex_particle {
			double X, Z;	// Global cartesian coordinate
			std::complex<double> x; // Foil eta plane coordinate
			double vX, vZ;	// Global velocities.
			double vorticity;	// vorticity of the particle.
		};
		std::vector<vortex_particle> m_vortex_particles;


	private:
		double particle_arg_eta_plain(const vortex_particle & particle);
		void compute_eta_plane_velocities(double time);
		void update_eta_plane_vortex_particle_coordinates(double time);
	};
} // End mFlow

