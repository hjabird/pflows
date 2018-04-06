#include "PlanarWakeULLT.h"


#include <HBTK/CartesianFiniteLine.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Checks.h>
#include <HBTK/Constants.h>
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


void mFlow::PlanarWakeULLT::wake_to_vtk(std::ostream &out_stream)
{
	auto wake = generate_planar_wake_object();
	wake.save_to_vtk(out_stream);
}

double mFlow::PlanarWakeULLT::compute_lift_coefficient()
{
	double coeff = 0;
	auto segments = segment_span_by_inner_solution();
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		double chord = wing_projection.chord(inner_solution_planes[i].origin().y());
		double section_width = segments[i].vector().magnitude();
		double cl_inner, cd_inner;
		std::tie(cl_inner, cd_inner) = inner_solutions[i].aerofoil_lift_and_drag_coefficients();
		coeff += cl_inner * section_width * chord;
	}
	coeff /= wing_projection.area();
	assert(HBTK::check_finite(coeff));
	return coeff;
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
	// Wing geometry (well, the lifting line's geometry):
	for (int i = 0; i < (int)y_positions.size(); i++) {
		wake.vertex(0, i, 
			HBTK::CartesianPoint2D({ 0.0, y_positions[i] })
		);
	}
	// And the vorticity of the wing
	for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
		wake.ring_strength(0, iy, inner_solutions[iy].bound_vorticity());
	}

	std::vector<double> wake_vorticity_acc(inner_solution_y_positions.size());
	for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
		wake_vorticity_acc[iy] = inner_solutions[iy].bound_vorticity();
	}
	// Wake geometry and vorticity:
	int num_wake_points = num_vortex_particles_per_inner_solution();
	std::vector<double> vortex_x_positions(inner_solution_y_positions.size());
	for (int ix = num_wake_points - 1; ix >= 0 ; ix--) {
		for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
			vortex_x_positions[iy] = inner_solutions[iy].m_vortex_particles[ix].position.x()
				- wing_projection.trailing_edge_X(inner_solution_planes[iy].origin().y());
		}
		// We use cubic spline interpolation. Edges might be dodgy - perhaps constrain
		// to convected position (Nope: otherwise vortex filament is unlikely to pass through vortex...)
		HBTK::CubicSpline1D line_spline(
			inner_solution_y_positions,
			vortex_x_positions);
		for (int iy = 0; iy < (int)y_positions.size(); iy++) {
			wake.vertex(num_wake_points - ix, iy,
				HBTK::CartesianPoint2D({ line_spline(y_positions[iy]), y_positions[iy] }));
		}
		if (ix > 0) {	// Fewer ring strengh values to set than vertices.
			for (int iy = 0; iy < (int)inner_solutions.size(); iy++) {
				// And use the spine's curvature to correct for the vorticity of the vortex ring.
				double angle = atan(line_spline.derivative(segments[iy].midpoint().y()));
				wake_vorticity_acc[iy] += inner_solutions[iy].m_vortex_particles[ix].vorticity
					/ cos(angle);
				wake.ring_strength(num_wake_points - ix, iy, wake_vorticity_acc[iy]);
			}
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
		double trailing_edge_x = wing_projection.trailing_edge_X(inner_solution_planes[i].origin().y());
		double downwash = 0;
		if (wake.size() > 0) {
			downwash += wake.patch_x_filament_downwash_inclusive(in_plane_coordinates[i],
				0, wake_depth, 0, wake_width);
			if (!quasi_steady) {
				downwash += wake.patch_y_filament_downwash_inclusive(in_plane_coordinates[i],
					1, wake_depth, 0, wake_width);
				// Now remove the W_wi as if were of infinite span (avoid W_wi in inner AND outer solution)
				for (int j = 0; j < (int)inner_solutions[i].m_vortex_particles.size(); j++) {
					auto & particle = inner_solutions[i].m_vortex_particles[j];
					downwash += particle.vorticity / (2 * HBTK::Constants::pi() * (particle.position.x() - trailing_edge_x));
				}
			}
		}
		if (!HBTK::check_finite(downwash)) {
			throw std::domain_error(
				"mFlow::PlanarWakeULLT::set_inner_solution_downwash: "
				"Computed non-finite downwash. "
				__FILE__ " : " + std::to_string(__LINE__)
			);
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
