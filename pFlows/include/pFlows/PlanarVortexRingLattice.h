#pragma once

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
		double patch_downwash_inclusive(const HBTK::CartesianPoint2D & measurement_position);
		// Compute downwash due to a patch in the lattice. The rings adjacent to 
		// patch influence the vorticity of the patch edge filaments.
		double patch_downwash_inclusive(const HBTK::CartesianPoint2D & measurement_position,
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