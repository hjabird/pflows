#pragma once
/*////////////////////////////////////////////////////////////////////////////
PlanarWakeULLT.h

An unsteady lifting theory based on Ramesh's LAUTAT in the inner solution.
Fails as amplitude inceases due to geometric non-linearity.

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

		// Vortex ring warping corrections - makes code blow up right now
		bool vortex_ring_warping_correction;
		// Symmetry
		bool symmetric;
		// The symmetry plane. 
		HBTK::CartesianPlane symmetry_plane;

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

		// The vortex ring strengths 
		std::vector<std::vector<double>> m_original_ring_strengths;

		// Correct for vortex filament curvature.
		void add_new_ring_to_ring_strengths();
		void update_inner_solution_vorticities_for_vortex_filament_curvature(PlanarVortexRingLattice & wake);

		// Turn the span of the wing into a set of segments, adding in symmetry if needed.
		std::vector<HBTK::CartesianFiniteLine2D> segment_span_by_inner_solution();

		// Subfunctions in wake vortex lattice generation
		std::vector<double> segment_endpoint_y_positions(const std::vector<HBTK::CartesianFiniteLine2D> & segments);
		std::vector<double> inner_solution_y_positions();
		void apply_lifting_line_vorticity(PlanarVortexRingLattice & lattice);
		void apply_lifting_line_geometry(PlanarVortexRingLattice & lattice, const std::vector<double> & y_positions);

		// The y_ordering of the inner solutions. The vector contains the index of an inner solution.
		// if the index is more than the number of solutions it is idx - num(inner_solutions) index
		// where the solution is reflected across the symmetry plane.
		std::vector<int> m_inner_solution_ordering;
		void recalculate_inner_solution_order();
		int reindexed_inner_solution(int idx) const;
		HBTK::CartesianPoint3D origin(int idx) const;
		double bound_vorticity(int idx);

		// Validity checks.
		bool valid_inner_solution_plane_count();
		bool matching_inner_solution_vortex_count();
		bool correct_ordering_of_inner_solution_planes();
	};
}