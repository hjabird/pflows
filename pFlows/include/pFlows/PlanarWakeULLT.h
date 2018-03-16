#pragma once


#include <functional>
#include <vector>

#include <HBTK/CartesianPlane.h>

#include "PlanarVortexRingLattice.h"
#include "Ramesh2014.h"
#include "VortexGroup2D.h"
#include "WingProjectionGeometry.h"

namespace mFlow {
	class LUALLT_planar_wake {
	public:
		LUALLT_planar_wake();
		~LUALLT_planar_wake();

		// Kinematics:
		std::function<double(double)> wing_Z;
		std::function<double(double)> wing_dZdt;
		std::function<double(double)> wing_AoA;
		std::function<double(double)> wing_dAoAdt;
		double pitch_location;
		// Time control:
		double delta_t;
		double time;

		// Wing projection on plane:
		mFlow::WingProjectionGeometry wing_projection;

		// Wake plane:
		HBTK::CartesianPlane wake_plane() const;

		// Inner solution planes:
		int number_of_inner_solution_planes;
		std::vector<HBTK::CartesianPlane> inner_solution_planes();
		std::vector<mFlow::VortexGroup2D> inner_solution_vortices();

		// Continue one time step
		void advance_one_step();

		// Compute coefficients for entire wing:
		double compute_lift_coefficient();
		double compute_drag_coefficient();
		double compute_moment_coeffient();

	private:
		// Inner solutions.
		// We expect them to be ordered with increasing y position.
		std::vector<Ramesh2014> m_inner_solutions;
		std::vector<HBTK::CartesianPlane> m_inner_solution_planes;

		// Wake
		PlanarVortexRingLattice generate_planar_wake_object();

		// Compute the downwash on each inner solution due to the wake and 
		// set the inner solution downwash variable to reflect this.
		void set_inner_solution_downwash(PlanarVortexRingLattice & wake);

		// Get the number of vortex particles in the inner solution.
		int num_vortex_particles_per_inner_solution();

		// Turn the span of the wing into a set of segments.
		std::vector<HBTK::CartesianFiniteLine2D> segment_span_by_inner_solution();

		// Validity checks.
		bool valid_inner_solution_plane_count();
		bool matching_inner_solution_vortex_count();
		bool correct_ordering_of_inner_solution_planes();
	};
}