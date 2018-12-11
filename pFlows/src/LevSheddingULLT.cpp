#include "LevSheddingULLT.h"
/*////////////////////////////////////////////////////////////////////////////
LevSheddingULLT.cpp

An unsteady lifting line theory utilisting a three dimensional wake model
with leading edge vortex shedding to inlude goemetric and aerodynamic 
non-linearity within a the lifting-line theory framework. Unsurprisingly,
this method does not work well.

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

#include <algorithm>
#include <array>
#include <numeric>

#include <HBTK/CartesianFiniteLine.h>
#include <HBTK/CartesianPoint.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Checks.h>
#include <HBTK/Constants.h>
#include <HBTK/CubicSpline1D.h>
#include <HBTK/CubicSplineND.h>
#include <HBTK/PotentialFlowDistributions.h>
#include <HBTK/StructuredBlockIndexerND.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>


mFlow::LevSheddingULLT::LevSheddingULLT()
	: quasi_steady(false),
	symmetric(false),
	vortex_ring_warping_correction(true)
{
}

mFlow::LevSheddingULLT::~LevSheddingULLT()
{
}

void mFlow::LevSheddingULLT::advance_one_step()
{
	VortexRingLattice tevwake = generate_tev_wake_object();
	VortexRingLattice levwake = generate_lev_wake_object();
	set_inner_solution_downwash(tevwake, levwake);
#pragma omp parallel for
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		inner_solutions[i].advance_one_step();
	}
}

void mFlow::LevSheddingULLT::initialise()
{
	for (auto & inner_sol : inner_solutions) {
		inner_sol.initialise();
		assert(inner_sol.m_te_vortex_particles.size() == 0);
	}
	for (auto & plane : inner_solution_planes) {
		double y_p = plane.origin().y();
		if (y_p >= wing_projection.semispan() || y_p <= -wing_projection.semispan()) {
			throw std::invalid_argument("In LevSheddingULLT::initialise: "
				"Inner solution wake plane is on wing tip or outside span. Not allowed. "
				__FILE__ " : " + std::to_string(__LINE__));
		}
	}

	recalculate_inner_solution_order();
}


void mFlow::LevSheddingULLT::wake_to_vtk(std::ostream &out_stream)
{
	auto twake = generate_tev_wake_object();
	auto lwake = generate_lev_wake_object();
	HBTK::Vtk::VtkUnstructuredDataset twd = twake.to_vtk_data();
	HBTK::Vtk::VtkUnstructuredDataset lwd = lwake.to_vtk_data();
	HBTK::Vtk::VtkWriter writer;
	writer.open_file(out_stream, HBTK::Vtk::VtkWriter::vtk_file_type::UnstructuredGrid);
	writer.appended = false;
	writer.write_piece(out_stream, twd);
	writer.write_piece(out_stream, lwd);
	writer.close_file(out_stream);
}

double mFlow::LevSheddingULLT::compute_lift_coefficient()
{
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

mFlow::VortexRingLattice mFlow::LevSheddingULLT::generate_tev_wake_object()
{
	VortexRingLattice wake(num_vortex_particles_per_inner_solution(),
		(int)m_inner_solution_ordering.size());
	if (wake.size() == 0) {
		return wake;
	}

	std::vector<HBTK::CartesianFiniteLine3D> segments = segment_span_by_inner_solution();
	std::vector<double> segment_end_y_positions = segment_endpoint_y_positions(segments);
	std::vector<double> inner_solution_y_positions = LevSheddingULLT::inner_solution_y_positions();

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
	std::vector<double> vortex_z_positions(inner_solution_y_positions.size());
	for (int ix = num_wake_points - 1; ix >= 0; ix--) {
		for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
			vortex_x_positions[iy] = inner_solutions[reindexed_inner_solution(iy)].m_te_vortex_particles[ix].position.x()
				// - inner_solutions[reindexed_inner_solution(iy)].foil_coordinate(1.).x();
				- wing_projection.trailing_edge_X(origin(iy).y());
			vortex_z_positions[iy] = inner_solutions[reindexed_inner_solution(iy)].m_te_vortex_particles[ix].position.y()
				- inner_solutions[reindexed_inner_solution(iy)].foil_coordinate(1.).y();
			// - wing_projection.trailing_edge_X(origin(iy).y());
		}
		// We use cubic spline interpolation. Edges might be dodgy - perhaps constrain
		// to convected position (Nope: otherwise vortex filament is unlikely to pass through vortex...)
		HBTK::CubicSpline1D line_spline_x(inner_solution_y_positions, vortex_x_positions);
		HBTK::CubicSpline1D line_spline_y(inner_solution_y_positions, inner_solution_y_positions);
		HBTK::CubicSpline1D line_spline_z(inner_solution_y_positions, vortex_z_positions);
		HBTK::CubicSplineND<3> spline3d({ line_spline_x, line_spline_y, line_spline_z });
		for (int iy = 0; iy < (int)segment_end_y_positions.size(); iy++) {
			wake.vertex(num_wake_points - ix, iy,
				HBTK::CartesianPoint3D(spline3d(segment_end_y_positions[iy])));
		}
		if (ix > 0) {	// Fewer ring strength values to set than vertices.
			for (int iy = 0; iy < (int)inner_solution_y_positions.size(); iy++) {
				wake_vorticity_acc[iy] += inner_solutions[reindexed_inner_solution(iy)].m_te_vortex_particles[ix].vorticity;
				wake.ring_strength(num_wake_points - ix, iy, wake_vorticity_acc[iy]);
			}
		}
	}
	return wake;
}

mFlow::VortexRingLattice mFlow::LevSheddingULLT::generate_lev_wake_object()
{
	VortexRingLattice wake(num_vortex_particles_per_inner_solution(),
		(int)m_inner_solution_ordering.size());
	if (wake.size() == 0) {
		return wake;
	}

	std::vector<HBTK::CartesianFiniteLine3D> segments = segment_span_by_inner_solution();
	std::vector<double> segment_end_y_positions = segment_endpoint_y_positions(segments);
	std::vector<double> inner_solution_y_positions = LevSheddingULLT::inner_solution_y_positions();

	// Lifting line geometry:
	apply_lifting_line_geometry(wake, segment_end_y_positions);
	apply_lifting_line_zero_vorticity(wake);

	std::vector<double> wake_vorticity_acc(inner_solution_y_positions.size());
	for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
		wake_vorticity_acc[iy] = 0.0;
	}
	// Wake geometry and vorticity:
	int num_wake_points = num_vortex_particles_per_inner_solution();
	std::vector<double> vortex_x_positions(inner_solution_y_positions.size());
	std::vector<double> vortex_z_positions(inner_solution_y_positions.size());
	for (int ix = num_wake_points - 1; ix >= 0; ix--) {
		for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
			vortex_x_positions[iy] = inner_solutions[reindexed_inner_solution(iy)].m_le_vortex_particles[ix].position.x()
				// - inner_solutions[reindexed_inner_solution(iy)].foil_coordinate(1.).x();
				- wing_projection.leading_edge_X(origin(iy).y());
			vortex_z_positions[iy] = inner_solutions[reindexed_inner_solution(iy)].m_le_vortex_particles[ix].position.y()
				- inner_solutions[reindexed_inner_solution(iy)].foil_coordinate(-1.).y();
			// - wing_projection.trailing_edge_X(origin(iy).y());
		}
		// We use cubic spline interpolation. Edges might be dodgy - perhaps constrain
		// to convected position (Nope: otherwise vortex filament is unlikely to pass through vortex...)
		HBTK::CubicSpline1D line_spline_x(inner_solution_y_positions, vortex_x_positions);
		HBTK::CubicSpline1D line_spline_y(inner_solution_y_positions, inner_solution_y_positions);
		HBTK::CubicSpline1D line_spline_z(inner_solution_y_positions, vortex_z_positions);
		HBTK::CubicSplineND<3> spline3d({ line_spline_x, line_spline_y, line_spline_z });
		for (int iy = 0; iy < (int)segment_end_y_positions.size(); iy++) {
			wake.vertex(num_wake_points - ix, iy,
				HBTK::CartesianPoint3D(spline3d(segment_end_y_positions[iy])));
		}
		if (ix > 0) {	// Fewer ring strength values to set than vertices.
			for (int iy = 0; iy < (int)inner_solution_y_positions.size(); iy++) {
				wake_vorticity_acc[iy] += inner_solutions[reindexed_inner_solution(iy)].m_le_vortex_particles[ix].vorticity;
				wake.ring_strength(num_wake_points - ix, iy, wake_vorticity_acc[iy]);
			}
		}
	}
	return wake;
}

void mFlow::LevSheddingULLT::set_inner_solution_downwash(VortexRingLattice & tev_wake, VortexRingLattice & lev_wake)
{
	std::vector<HBTK::CartesianPoint3D> coordinates;
	for (auto & plane : inner_solution_planes) {
		coordinates.push_back(plane.origin());
	}
	int wake_depth = tev_wake.extent()[0];
	int wake_width = tev_wake.extent()[1];
	for (int i = 0; i < (int)inner_solutions.size(); i++) {
		double trailing_edge_x = wing_projection.trailing_edge_X(inner_solution_planes[i].origin().y());
		HBTK::CartesianVector3D downwash({ 0, 0, 0 });
		if (tev_wake.size() > 0) {
			downwash += tev_wake.patch_x_filament_induced_velocity_inclusive(coordinates[i],
				0, wake_depth, 0, wake_width);
			downwash += lev_wake.patch_x_filament_induced_velocity_inclusive(coordinates[i],
				0, wake_depth, 0, wake_width);
			if (!quasi_steady) {
				downwash += tev_wake.patch_y_filament_induced_velocity_inclusive(coordinates[i],
					1, wake_depth, 0, wake_width);
				downwash += lev_wake.patch_y_filament_induced_velocity_inclusive(coordinates[i],
					1, wake_depth, 0, wake_width);
				// Now remove the W_wi as if were of infinite span (avoid W_wi in inner AND outer solution)
				for (int j = 0; j < (int)inner_solutions[i].m_te_vortex_particles.size(); j++) {
					auto & tev_particle = inner_solutions[i].m_te_vortex_particles[j];
					HBTK::CartesianVector2D diff = tev_particle.position - inner_solutions[i].foil_coordinate(1);
					downwash.x() += tev_particle.vorticity * HBTK::PointVortex::unity_u_vel(0.0, 0.0, diff.x(), diff.y());
					downwash.z() += tev_particle.vorticity * HBTK::PointVortex::unity_v_vel(0.0, 0.0, diff.x(), diff.y());
					auto & lev_particle = inner_solutions[i].m_le_vortex_particles[j];
					diff = lev_particle.position - inner_solutions[i].foil_coordinate(-1);
					downwash.x() += lev_particle.vorticity * HBTK::PointVortex::unity_u_vel(0.0, 0.0, diff.x(), diff.y());
					downwash.z() += lev_particle.vorticity * HBTK::PointVortex::unity_v_vel(0.0, 0.0, diff.x(), diff.y());
				}
			}
		}
		if (!HBTK::check_finite(downwash)) {
			throw std::domain_error(
				"mFlow::LevSheddingULLT::set_inner_solution_downwash: "
				"Computed non-finite downwash. "
				__FILE__ " : " + std::to_string(__LINE__)
			);
		}
		inner_solutions[i].external_purturbation =
			[=](HBTK::CartesianPoint2D p, double t) { return HBTK::CartesianVector2D({ downwash.x(), downwash.z() }); };
	}
	return;
}

int mFlow::LevSheddingULLT::num_vortex_particles_per_inner_solution()
{
	return (int)inner_solutions[0].m_te_vortex_particles.size();
}

std::vector<HBTK::CartesianFiniteLine3D> mFlow::LevSheddingULLT::segment_span_by_inner_solution()
{
	std::vector<HBTK::CartesianFiniteLine3D> segments((int)m_inner_solution_ordering.size());
	segments[0].start() = HBTK::CartesianPoint3D({ 0, -wing_projection.semispan(), 0 });
	for (int i = 1; i < (int)m_inner_solution_ordering.size(); i++) {
		double y_pos = (origin(i - 1).y() + origin(i).y()) / 2;	// Average
		segments[i - 1].end() = HBTK::CartesianPoint3D({ 0, y_pos, 0 });
		segments[i].start() = HBTK::CartesianPoint3D({ 0, y_pos, 0 });
	}
	segments.back().end() = HBTK::CartesianPoint3D({ 0, wing_projection.semispan(), 0. });
	assert(HBTK::check_finite(segments));
	return segments;
}

std::vector<double> mFlow::LevSheddingULLT::segment_endpoint_y_positions(const std::vector<HBTK::CartesianFiniteLine3D>& segments)
{
	std::vector<double> y_positions;
	for (auto & segment : segments) {
		y_positions.emplace_back(segment.start().y());
	}
	y_positions.emplace_back(segments.back().end().y());
	return y_positions;
}

std::vector<double> mFlow::LevSheddingULLT::inner_solution_y_positions()
{
	std::vector<double> y;
	for (int i = 0; i < (int)m_inner_solution_ordering.size(); i++) {
		y.push_back(origin(i).y());
	}
	return y;
}

void mFlow::LevSheddingULLT::apply_lifting_line_vorticity(VortexRingLattice & lattice)
{
	for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
		lattice.ring_strength(0, iy, bound_vorticity(iy));
	}
	return;
}

void mFlow::LevSheddingULLT::apply_lifting_line_geometry(VortexRingLattice & lattice, const std::vector<double> & y_positions)
{
	std::vector<double> inner_solution_y_positions;
	// Wing geometry (well, the lifting line's geometry):
	for (int i = 0; i < (int)y_positions.size(); i++) {
		lattice.vertex(0, i,
			HBTK::CartesianPoint3D({ 0.0, y_positions[i], 0.0 })
		);
	}
}

void mFlow::LevSheddingULLT::apply_lifting_line_zero_vorticity(VortexRingLattice & lattice)
{
	for (int iy = 0; iy < (int)m_inner_solution_ordering.size(); iy++) {
		lattice.ring_strength(0, iy, 0.0);
	}
	return;
}

void mFlow::LevSheddingULLT::recalculate_inner_solution_order()
{
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

int mFlow::LevSheddingULLT::reindexed_inner_solution(int idx) const
{
	assert(idx >= 0);
	assert(idx < (int)m_inner_solution_ordering.size());
	return m_inner_solution_ordering[idx] % inner_solutions.size();
}

HBTK::CartesianPoint3D mFlow::LevSheddingULLT::origin(int idx) const
{
	HBTK::CartesianPoint3D pnt;
	pnt = inner_solution_planes[reindexed_inner_solution(idx)].origin();
	if (idx >= (int)inner_solutions.size()) {
		pnt = symmetry_plane.symmetric_point(pnt);
	}
	return pnt;
}

double mFlow::LevSheddingULLT::bound_vorticity(int idx)
{
	return inner_solutions[reindexed_inner_solution(idx)].bound_vorticity();
}

bool mFlow::LevSheddingULLT::valid_inner_solution_plane_count()
{
	return ((int)inner_solutions.size() > 0 ? true : false);
}

bool mFlow::LevSheddingULLT::matching_inner_solution_vortex_count()
{
	int count = inner_solutions[0].m_te_vortex_particles.size();
	for (auto & solution : inner_solutions) {
		if (count != solution.m_te_vortex_particles.size()) return false;
	}
	return true;
}

bool mFlow::LevSheddingULLT::correct_ordering_of_inner_solution_planes()
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

