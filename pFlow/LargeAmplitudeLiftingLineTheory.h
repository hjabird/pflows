#pragma once

#include <vector>

#include <HBTK\CartesianPoint.h>
#include <HBTK\CartesianPlane.h>

#include "WingProjectionGeometry.h"

namespace mFlow {
	class LargeAmplitudeLiftingLineTheory
	{
	public:
		LargeAmplitudeLiftingLineTheory();
		~LargeAmplitudeLiftingLineTheory();

		WingProjectionGeometry wing;

		struct vortex {
			HBTK::CartesianPoint2D location;
			double vorticity;
		};
	private:

		// Our vortex particles are defined by their 2D positions in planes:
		std::vector<HBTK::CartesianPlane> m_planes;		// (Assumed ordered -y to y)
		// With the particles defined as:
		std::vector<std::vector<vortex>> m_vortex_particles;
		// Initialise the position of the planes:
		void initialise_inner_solution_planes(std::vector<double> y_locations);

		// The chord lengths at the positions on the plane:
		std::vector<double> m_chord_lengths;
		// Initialised by (once we have initialised plane positions and our wing)
		void initialise_chord_lengths();

		std::vector<double> get_upwash();
		std::vector<HBTK::CartesianLine3D> get_segments();

		bool check_matching_vortex_particles();
	};
}
