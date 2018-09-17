#pragma once
/*////////////////////////////////////////////////////////////////////////////
VortexRingLattice.cpp

A singular vortex lattice representation.

Copyright 2017-2018 HJA Bird

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
#include <HBTK/VtkUnstructuredDataset.h>

namespace mFlow {
	class VortexRingLattice {
	public:
		VortexRingLattice(int x_extent, int y_extent);
		~VortexRingLattice();

		// Compute the downwash at a point in the plane due to the 
		// entire lattice.
		HBTK::CartesianVector3D induced_velocity(const HBTK::CartesianPoint3D & measurement_position);
		// Compute downwash due to a patch in the lattice. The rings adjacent to 
		// patch influence the vorticity of the patch edge filaments.
		HBTK::CartesianVector3D patch_induced_velocity_inclusive(const HBTK::CartesianPoint3D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);
		HBTK::CartesianVector3D patch_x_filament_induced_velocity_inclusive(const HBTK::CartesianPoint3D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);
		HBTK::CartesianVector3D patch_y_filament_induced_velocity_inclusive(const HBTK::CartesianPoint3D & measurement_position,
			int min_x, int max_x, int min_y, int max_y);

		// Methods to extract vortex filaments:
		// For index in [0, nx), [0, ny] 
		HBTK::CartesianFiniteLine3D edge_x(int index_x, int index_y);
		// For index in [0, nx), [0, ny]
		double edge_x_vorticity(int index_x, int index_y);
		// For index in [0, nx], [0, ny)
		HBTK::CartesianFiniteLine3D edge_y(int index_x, int index_y);
		// For index in [0, nx], [0, ny)
		double edge_y_vorticity(int index_x, int index_y);

		// Methods to set geometry and ring strength:
		// For index in [0, nx), [0, ny)
		double ring_strength(int index_x, int index_y);
		void ring_strength(int index_x, int index_y, double vorticity);
		HBTK::CartesianPoint3D vertex(int index_x, int index_y);
		void vertex(int index_x, int index_y, HBTK::CartesianPoint3D position);

		std::array<int, 2> extent() const;
		int size() const;

		HBTK::Vtk::VtkUnstructuredDataset to_vtk_data();

	private:
		int m_extent_x;
		int m_extent_y;

		// nx+1 by ny+1 structured block of the ring corners
		HBTK::StructuredValueBlockND<2, HBTK::CartesianPoint3D> m_geometry;
		// nx by ny structured block of the ring vorticities.
		HBTK::StructuredValueBlockND<2, double> m_ring_strenghts;

		// Down wash from a vortex filament (all in same plane)
		HBTK::CartesianVector3D filament_induced_velocity(
			const HBTK::CartesianFiniteLine3D & filament,
			double strength,
			const HBTK::CartesianPoint3D & measurement_point
		);

	};
}