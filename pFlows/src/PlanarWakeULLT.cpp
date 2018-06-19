#include "PlanarWakeULLT.h"

#include <algorithm>
#include <numeric>

#include <HBTK/CartesianFiniteLine.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Checks.h>
#include <HBTK/Constants.h>
#include <HBTK/CubicSpline1D.h>
#include <HBTK/StructuredBlockIndexerND.h>
#include <HBTK/RuntimeProfiler.h>


mFlow::PlanarWakeULLT::PlanarWakeULLT()
	: quasi_steady(false),
	symmetric(false),
	vortex_ring_warping_correction(true)
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
	HBTK::RuntimeProfiler(__FUNCTION__, __LINE__, true);
	PlanarVortexRingLattice wake = generate_planar_wake_object();
	update_inner_solution_vorticities_for_vortex_filament_curvature(wake);
	set_inner_solution_downwash(wake);
#pragma omp parallel for
	for (int i = 0; i < (int) inner_solutions.size(); i++) {
		inner_solutions[i].advance_one_step();
	}
	add_new_ring_to_ring_strengths();
}

void mFlow::PlanarWakeULLT::initialise()
{
	for (auto & inner_sol : inner_solutions) {
		inner_sol.initialise();
		assert(inner_sol.m_te_vortex_particles.size() == 0);
	}
	recalculate_inner_solution_order();
	m_original_ring_strengths.resize(inner_solutions.size());
}


void mFlow::PlanarWakeULLT::wake_to_vtk(std::ostream &out_stream)
{
	auto wake = generate_planar_wake_object();
	wake.save_to_vtk(out_stream);
}

double mFlow::PlanarWakeULLT::compute_lift_coefficient()
{
	HBTK::RuntimeProfiler(__FUNCTION__, __LINE__, true);
	double coeff = 0;
	auto segments = segment_span_by_inner_solution();
	std::vector<double> inner_solution_ys = inner_solution_y_positions();
#pragma omp parallel for
	for (int i = 0; i < (int)m_inner_solution_ordering.size(); i++) {
		double chord = wing_projection.chord(inner_solution_ys[i]);
		double section_width = segments[i].vector().magnitude();
		double cl_inner, cd_inner;
		std::tie(cl_inner, cd_inner) = inner_solutions[reindexed_inner_solution(i)].aerofoil_lift_and_drag_coefficients();
#pragma omp critical
		coeff += cl_inner * section_width * chord;
	}
	coeff /= wing_projection.area();
	assert(HBTK::check_finite(coeff));
	return coeff;
}

mFlow::PlanarVortexRingLattice mFlow::PlanarWakeULLT::generate_planar_wake_object()
{
	HBTK::RuntimeProfiler(__FUNCTION__, __LINE__, true);
	PlanarVortexRingLattice wake(num_vortex_particles_per_inner_solution(),
		(int)m_inner_solution_ordering.size());
	if (wake.size() == 0) {
		return wake;
	}

	std::vector<HBTK::CartesianFiniteLine2D> segments = segment_span_by_inner_solution();
	std::vector<double> segment_end_y_positions	= segment_endpoint_y_positions(segments);
	std::vector<double> inner_solution_y_positions = PlanarWakeULLT::inner_solution_y_positions();

	// Lifting line geometry:
	apply_lifting_line_geometry(wake, segment_end_y_positions);
	apply_lifting_line_vorticity(wake);

	std::vector<double> wake_vorticity_acc(inner_solution_y_positions.size());
	for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
		wake_vorticity_acc[iy] = bound_vorticity(iy);
	}
	// Wake geometry and vorticity:
	int num_wake_points = num_vortex_particles_per_inner_solution();
	std::vector<double> vortex_x_positions(inner_solution_y_positions.size());
	for (int ix = num_wake_points - 1; ix >= 0 ; ix--) {
		for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
			vortex_x_positions[iy] = inner_solutions[reindexed_inner_solution(iy)].m_te_vortex_particles[ix].position.x()
				// - inner_solutions[reindexed_inner_solution(iy)].foil_coordinate(1.).x();
				- wing_projection.trailing_edge_X(origin(iy).y());
		}
		// We use cubic spline interpolation. Edges might be dodgy - perhaps constrain
		// to convected position (Nope: otherwise vortex filament is unlikely to pass through vortex...)
		HBTK::CubicSpline1D line_spline(
			inner_solution_y_positions,
			vortex_x_positions);
		for (int iy = 0; iy < (int)segment_end_y_positions.size(); iy++) {
			wake.vertex(num_wake_points - ix, iy,
				HBTK::CartesianPoint2D({ line_spline(segment_end_y_positions[iy]), segment_end_y_positions[iy] }));
		}
		if (ix > 0) {	// Fewer ring strength values to set than vertices.
			for (int iy = 0; iy < (int)inner_solution_y_positions.size(); iy++) {
				// And use the spine's curvature to correct for the vorticity of the vortex ring.
				double angle = atan(line_spline.derivative(segments[iy].midpoint().y()));
				wake_vorticity_acc[iy] += m_original_ring_strengths[reindexed_inner_solution(iy)][ix];
				wake.ring_strength(num_wake_points - ix, iy, wake_vorticity_acc[iy]);
			}
		}
	}
	return wake;
}

void mFlow::PlanarWakeULLT::set_inner_solution_downwash(PlanarVortexRingLattice & wake)
{
	HBTK::RuntimeProfiler(__FUNCTION__, __LINE__, true);
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
				for (int j = 0; j < (int)inner_solutions[i].m_te_vortex_particles.size(); j++) {
					auto & particle = inner_solutions[i].m_te_vortex_particles[j];
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
	return (int) inner_solutions[0].m_te_vortex_particles.size();
}

std::vector<HBTK::CartesianFiniteLine2D> mFlow::PlanarWakeULLT::segment_span_by_inner_solution()
{
	// assert(correct_ordering_of_inner_solution_planes());

	std::vector<HBTK::CartesianFiniteLine2D> segments((int) m_inner_solution_ordering.size());
	segments[0].start() = HBTK::CartesianPoint2D({ 0, -wing_projection.semispan() });
	for (int i = 1; i < (int)m_inner_solution_ordering.size(); i++) {
		double y_pos = (origin(i - 1).y() + origin(i).y()) / 2;	// Average
		segments[i - 1].end() = HBTK::CartesianPoint2D({ 0, y_pos });
		segments[i].start() = HBTK::CartesianPoint2D({ 0, y_pos });
	}
	segments.back().end() = HBTK::CartesianPoint2D({ 0, wing_projection.semispan() });

	assert(HBTK::check_finite(segments));
	return segments;
}

std::vector<double> mFlow::PlanarWakeULLT::segment_endpoint_y_positions(const std::vector<HBTK::CartesianFiniteLine2D>& segments)
{
	std::vector<double> y_positions;
	for (auto & segment : segments) {
		y_positions.emplace_back(segment.start().y());
	}
	y_positions.emplace_back(segments.back().end().y());
	return y_positions;
}

std::vector<double> mFlow::PlanarWakeULLT::inner_solution_y_positions()
{
	std::vector<double> y;
	for (int i = 0; i < (int)m_inner_solution_ordering.size(); i++) {
		y.push_back(origin(i).y());
	}
	return y;
}

void mFlow::PlanarWakeULLT::apply_lifting_line_vorticity(PlanarVortexRingLattice & lattice)
{
	for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
		lattice.ring_strength(0, iy, bound_vorticity(iy));
	}
	return;
}

void mFlow::PlanarWakeULLT::apply_lifting_line_geometry(PlanarVortexRingLattice & lattice, const std::vector<double> & y_positions)
{
	std::vector<double> inner_solution_y_positions;
	// Wing geometry (well, the lifting line's geometry):
	for (int i = 0; i < (int)y_positions.size(); i++) {
		lattice.vertex(0, i,
			HBTK::CartesianPoint2D({ 0.0, y_positions[i] })
		);
	}
}

void mFlow::PlanarWakeULLT::recalculate_inner_solution_order()
{
	HBTK::RuntimeProfiler( __FUNCTION__, __LINE__, true);
	std::vector<double> y_coords;
	for (auto & sol : inner_solution_planes) {
		y_coords.push_back(sol.origin().y());
	}
	if (symmetric) {
		for (auto & sol : inner_solution_planes) {
			y_coords.push_back(symmetry_plane.symmetric_point(sol.origin()).y());
		}
	}
	m_inner_solution_ordering.resize(y_coords.size());
	std::iota(m_inner_solution_ordering.begin(), m_inner_solution_ordering.end(), 0);
	// Sort the inner solution ordering by the y_coordinates.
	std::sort(m_inner_solution_ordering.begin(), m_inner_solution_ordering.end(),
		[&y_coords](int i1, int i2) { return y_coords[i1] < y_coords[i2]; });
	return;
}

int mFlow::PlanarWakeULLT::reindexed_inner_solution(int idx) const
{
	assert(idx >= 0);
	assert(idx < (int)m_inner_solution_ordering.size());
	return m_inner_solution_ordering[idx] % inner_solutions.size();
}

HBTK::CartesianPoint3D mFlow::PlanarWakeULLT::origin(int idx) const
{
	HBTK::CartesianPoint3D pnt;
	pnt = inner_solution_planes[reindexed_inner_solution(idx)].origin();
	if (idx >= (int)inner_solutions.size()) {
		pnt = symmetry_plane.symmetric_point(pnt);
	}
	return pnt;
}

double mFlow::PlanarWakeULLT::bound_vorticity(int idx)
{
	return inner_solutions[reindexed_inner_solution(idx)].bound_vorticity();
}

bool mFlow::PlanarWakeULLT::valid_inner_solution_plane_count()
{
	return ((int)inner_solutions.size() > 0 ? true : false);
}

bool mFlow::PlanarWakeULLT::matching_inner_solution_vortex_count()
{
	int count = inner_solutions[0].m_te_vortex_particles.size();
	for (auto & solution : inner_solutions) {
		if (count != solution.m_te_vortex_particles.size()) return false;
	}
	return true;
}

bool mFlow::PlanarWakeULLT::correct_ordering_of_inner_solution_planes()
{
	bool good;
	std::vector<double> inner_solution_y = inner_solution_y_positions();
	good = std::is_sorted(m_inner_solution_ordering.begin(), m_inner_solution_ordering.end(),
		[&inner_solution_y](int a, int b) {
		double y_a = inner_solution_y[a];
		double y_b = inner_solution_y[b];
		return y_a < y_b;
	});
	return good;
}


void mFlow::PlanarWakeULLT::add_new_ring_to_ring_strengths(void) 
{
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		m_original_ring_strengths[i].push_back(
			inner_solutions[i].m_te_vortex_particles.most_recently_added().vorticity);
	}
	return;
}

void mFlow::PlanarWakeULLT::update_inner_solution_vorticities_for_vortex_filament_curvature(PlanarVortexRingLattice & wake)
{
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		int y_idx = *(std::find(m_inner_solution_ordering.begin(),
			m_inner_solution_ordering.end(), i));
		for (int ix = 0; ix < num_vortex_particles_per_inner_solution(); ix++) {
			HBTK::CartesianVector2D tangent = wake.edge_y(ix + 1, y_idx).vector();
			if (!vortex_ring_warping_correction) {
				inner_solutions[i].m_te_vortex_particles[ix].vorticity =
					m_original_ring_strengths[i][ix];
			}
			else {
				inner_solutions[i].m_te_vortex_particles[ix].vorticity =
					m_original_ring_strengths[i][ix] * abs(tangent.y()) / tangent.magnitude();
			}
		}
	}
}
