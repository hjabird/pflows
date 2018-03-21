#pragma once


#include <functional>
#include <string>
#include <vector>

#include <HBTK/CartesianPlane.h>

#include "PlanarVortexRingLattice.h"
#include "Ramesh2014.h"
#include "VortexGroup2D.h"
#include "WingProjectionGeometry.h"

namespace mFlow {
	class PlanarWakeULLT {
	public:
		PlanarWakeULLT();
		~PlanarWakeULLT();

		// Kinematics and time control is defined within inner solutions.

		// Wing projection on plane:
		mFlow::WingProjectionGeometry wing_projection;

		// Wake plane:
		HBTK::CartesianPlane wake_plane() const;

		// Inner solution planes:
		std::vector<HBTK::CartesianPlane> inner_solution_planes;
		std::vector<mFlow::Ramesh2014> inner_solutions;

		// Continue one time step
		void advance_one_step();
		void initialise();

		// Who would want Quasi-steady? Default false.
		bool quasi_steady;

		// Save wake to file:
		void wake_to_vtk(std::ostream & out_stream);

		// Compute the lift coefficient
		double compute_lift_coefficient();

	private:
		// Wake:
		// Wing is on the x (idx0) = 0 line.
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