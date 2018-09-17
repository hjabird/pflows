#pragma once
/*////////////////////////////////////////////////////////////////////////////
PlanarVortexRingLattice.h

A vortex ring lattice constrained to a the x-y plane.

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

#include <string>

#include <HBTK/CartesianFiniteLine.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianRectilinearPanel.h>
#include <HBTK/StructuredValueBlockND.h>

namespace mFlow {
	class PlanarVortexRingLattice {
	public:
		PlanarVortexRingLattice(int x_extent, int y_extent);
		~PlanarVortexRingLattice();

		// Compute the downwash at a point in the plane due to the 
		// entire lattice.
		double downwash(const HBTK::CartesianPoint2D & measurement_position);
		// Compute downwash due to a patch in the lattice. The rings adjacent to 
		// patch influence the vorticity of the patch edge filaments.
		double patch_downwash_inclusive(const HBTK::CartesianPoint2D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);
		double patch_x_filament_downwash_inclusive(const HBTK::CartesianPoint2D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);
		double patch_y_filament_downwash_inclusive(const HBTK::CartesianPoint2D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);

		// Methods to extract vortex filaments:
		// For index in [0, nx), [0, ny] 
		HBTK::CartesianFiniteLine2D edge_x(int index_x, int index_y);
		// For index in [0, nx), [0, ny]
		double edge_x_vorticity(int index_x, int index_y);
		// For index in [0, nx], [0, ny)
		HBTK::CartesianFiniteLine2D edge_y(int index_x, int index_y);
		// For index in [0, nx], [0, ny)
		double edge_y_vorticity(int index_x, int index_y);

		// Methods to set geometry and ring strength:
		// For index in [0, nx), [0, ny)
		double ring_strength(int index_x, int index_y);
		void ring_strength(int index_x, int index_y, double vorticity);
		HBTK::CartesianPoint2D vertex(int index_x, int index_y);
		void vertex(int index_x, int index_y, HBTK::CartesianPoint2D position);
		
		std::array<int, 2> extent() const;
		int size() const;

		void save_to_vtk(std::ostream & out_stream);

	private:
		int m_extent_x;
		int m_extent_y;

		// nx+1 by ny+1 structured block of the ring corners
		HBTK::StructuredValueBlockND<2, HBTK::CartesianPoint2D> m_geometry;
		// nx by ny structured block of the ring vorticities.
		HBTK::StructuredValueBlockND<2, double> m_ring_strenghts;

		// Down wash from a vortex filament (all in same plane)
		double filament_downwash(
			const HBTK::CartesianFiniteLine2D & filament,
			double strength,
			const HBTK::CartesianPoint2D & measurement_point
		);

	};
}