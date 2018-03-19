#include "PlanarWakeULLT.h"


#include <HBTK/Checks.h>
#include <HBTK/CubicSpline1D.h>
#include <HBTK/StructuredBlockIndexerND.h>

mFlow::PlanarWakeULLT::PlanarWakeULLT()
	: quasi_steady(false)
{
}

mFlow::PlanarWakeULLT::~PlanarWakeULLT()
{
}

HBTK::CartesianPlane mFlow::PlanarWakeULLT::wake_plane() const
{
	HBTK::CartesianPlane plane(HBTK::CartesianPoint3D({0, 0, 0}), 
								HBTK::CartesianPoint3D({ 1, 0, 0 }), 
								HBTK::CartesianPoint3D({ 0, 1, 0 }));
	return plane;
}

void mFlow::PlanarWakeULLT::advance_one_step()
{
	PlanarVortexRingLattice wake = generate_planar_wake_object();
	set_inner_solution_downwash(wake);
#pragma omp parallel for
	for (int i = 0; i < (int) inner_solutions.size(); i++) {
		inner_solutions[i].advance_one_step();
	}
}

void mFlow::PlanarWakeULLT::initialise()
{
	for (auto & inner_sol : inner_solutions) {
		inner_sol.initialise();
	}
}


mFlow::PlanarVortexRingLattice mFlow::PlanarWakeULLT::generate_planar_wake_object()
{
	PlanarVortexRingLattice wake(num_vortex_particles_per_inner_solution(),
		(int)inner_solutions.size());
	if (wake.size() == 0) {
		return wake;
	}

	std::vector<HBTK::CartesianFiniteLine2D> segments = segment_span_by_inner_solution();
	std::vector<double> y_positions;
	for (auto & segment : segments) y_positions.emplace_back(segment.start().y());
	y_positions.emplace_back(segments.back().end().y());
	std::vector<double> inner_solution_y_positions;
	for (auto & plane : inner_solution_planes) inner_solution_y_positions.push_back(plane.origin().y());
	// Wing geometry:
	for (int i = 0; i < (int)y_positions.size(); i++) {
		wake.vertex(0, i, HBTK::CartesianPoint2D({ 0, y_positions[i] }));
	}
	// Wake geometry:
	int num_wake_points = num_vortex_particles_per_inner_solution();
	std::vector<double> vortex_x_positions(inner_solution_y_positions.size());
	for (int ix = num_wake_points - 1; ix >= 0 ; ix--) {
		for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
			vortex_x_positions[iy] = inner_solutions[iy].m_vortex_particles[ix].position.x()
				- inner_solution_planes[iy](
					HBTK::CartesianPoint2D({ wing_projection.semichord(y_positions[iy]), 0 })).x();
		}
		// We use cubic spline interpolation. Edges might be dodgy...
		HBTK::CubicSpline1D line_spline(
			inner_solution_y_positions,
			vortex_x_positions);
		for (int iy = 0; iy < (int)y_positions.size(); iy++) {
			wake.vertex(ix + 1, iy, 
				HBTK::CartesianPoint2D({ line_spline(y_positions[iy]), y_positions[iy] }));
		}
	}
	// Wake vorticity
	for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
		wake.ring_strength(0, iy, inner_solutions[iy].bound_vorticity());
	}
	for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
		double vorticity_acc = inner_solutions[iy].bound_vorticity();
		for (int ix = 0; ix < inner_solutions[iy].m_vortex_particles.size() - 1; ix++) {
			// Last vortex is correct due to Kelvin condition.
			vorticity_acc += inner_solutions[iy].m_vortex_particles[ix].vorticity;
			wake.ring_strength(ix + 1, iy, vorticity_acc);
		}
	}
	return wake;
}

void mFlow::PlanarWakeULLT::set_inner_solution_downwash(PlanarVortexRingLattice & wake)
{
	std::vector<HBTK::CartesianPoint3D> coordinates;
	for (auto & plane : inner_solution_planes) {
		coordinates.push_back(plane.origin());
	}
	std::vector<HBTK::CartesianPoint2D> in_plane_coordinates;
	for (auto & origin : coordinates) {
		in_plane_coordinates.push_back(wake_plane().projection(origin));
	}
	int wake_depth = wake.extent()[0];
	int wake_width = wake.extent()[1];
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		double downwash = 0;
		if (wake.size() > 0) {
			downwash += wake.patch_x_filament_downwash_inclusive(in_plane_coordinates[i],
				0, wake_depth, 0, wake_width);
			if (!quasi_steady) {
				downwash += wake.patch_y_filament_downwash_inclusive(in_plane_coordinates[i],
					1, wake_depth, 0, i);
				downwash += wake.patch_y_filament_downwash_inclusive(in_plane_coordinates[i],
					1, wake_depth, i + 1, wake_width);
			}
		}
		inner_solutions[i].free_stream_velocity.y() = downwash;
	}
	return;
}

int mFlow::PlanarWakeULLT::num_vortex_particles_per_inner_solution()
{
	return (int) inner_solutions[0].m_vortex_particles.size();
}

std::vector<HBTK::CartesianFiniteLine2D> mFlow::PlanarWakeULLT::segment_span_by_inner_solution()
{
	assert(correct_ordering_of_inner_solution_planes());
	std::vector<HBTK::CartesianFiniteLine2D> segments((int) inner_solutions.size());
	segments[0].start() = HBTK::CartesianPoint2D({ 0, -wing_projection.semispan() });
	for (int i = 1; i < (int)inner_solutions.size(); i++) {
		double y_pos = (inner_solution_planes[i - 1].origin().y() 
			+ inner_solution_planes[i].origin().y()) / 2;	// Average
		segments[i - 1].end() = HBTK::CartesianPoint2D({ 0, y_pos });
		segments[i].start() = HBTK::CartesianPoint2D({ 0, y_pos });
	}
	segments.back().end() = HBTK::CartesianPoint2D({ 0, wing_projection.semispan() });
	assert(HBTK::check_finite(segments));
	return segments;
}

bool mFlow::PlanarWakeULLT::valid_inner_solution_plane_count()
{
	return ((int)inner_solutions.size() > 0 ? true : false);
}

bool mFlow::PlanarWakeULLT::matching_inner_solution_vortex_count()
{
	int count = inner_solutions[0].m_vortex_particles.size();
	for (auto & solution : inner_solutions) {
		if (count != solution.m_vortex_particles.size()) return false;
	}
	return true;
}

bool mFlow::PlanarWakeULLT::correct_ordering_of_inner_solution_planes()
{
	double last_plane_y = inner_solution_planes[0].origin().y();
	bool good = true;
	for (int i = 1; i < (int)inner_solution_planes.size(); i++) {
		double y_pos = inner_solution_planes[i].origin().y();
		if (last_plane_y >= y_pos) good = false;
		last_plane_y = y_pos;
	}
	return good;
}
