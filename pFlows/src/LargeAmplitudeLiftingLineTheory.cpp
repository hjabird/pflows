#include "LargeAmplitudeLiftingLineTheory.h"

#include <cassert>

#include <HBTK/CartesianLine.h>
#include <HBTK/PotentialFlowDistributions.h>

mFlow::LargeAmplitudeLiftingLineTheory::LargeAmplitudeLiftingLineTheory()
{
}


mFlow::LargeAmplitudeLiftingLineTheory::~LargeAmplitudeLiftingLineTheory()
{
}

void mFlow::LargeAmplitudeLiftingLineTheory::initialise_inner_solution_planes(std::vector<double> y_locations)
{
	m_planes.resize(y_locations.size());
	for (int i = 0; i < (int) y_locations.size(); i++) {
		HBTK::CartesianVector3D offset({ 0, y_locations[i], 0 });
		HBTK::CartesianPoint3D origin({ 0,0,0 }), 
			primary_dir({ 1,0,0 }), secondary_dir({ 0,0,1 });
		origin = origin + offset;
		primary_dir = primary_dir + offset;
		secondary_dir = secondary_dir + offset;
		m_planes[i] = HBTK::CartesianPlane(origin, primary_dir, secondary_dir);
	}
}

void mFlow::LargeAmplitudeLiftingLineTheory::initialise_chord_lengths()
{
	m_chord_lengths.resize(m_planes.size());
	for (int i = 0; i < (int) m_planes.size(); i++) {
		HBTK::CartesianPoint3D plane_origin = m_planes[i](HBTK::CartesianPoint2D({ 0,0 }));
		double chord = wing.chord(plane_origin.y());
		m_chord_lengths[i] = chord;
	}
}

std::vector<double> mFlow::LargeAmplitudeLiftingLineTheory::get_upwash()
{
	assert(check_matching_vortex_particles());
	int n_planes = (int)m_planes.size();
	int n_particles_inner = (int)m_vortex_particles[0].size();
	std::vector<HBTK::CartesianVector3D> wash(n_planes);
	std::vector<HBTK::CartesianPoint3D> collocation_points(n_planes);
	for (int i = 0; i < n_planes; i++) { collocation_points[i] = m_planes[i].origin(); }
	std::vector<HBTK::CartesianLine3D> spanwise_segments = get_segments();

	// Consideration of the effect of the y dir part of the vortex rings.
	// Spanwise vortex ring consideration
	for (int i = 0; i < n_planes; i++) {
		HBTK::CartesianLine3D segment = spanwise_segments[i];
		// Wakewise vortex ring consideration
		for (int j = 0; j < n_particles_inner; j++) {
			double vorticity = m_vortex_particles[i][j].vorticity;
			HBTK::CartesianPoint3D vp = m_planes[i](m_vortex_particles[i][j].location);
			HBTK::CartesianVector3D displ = vp - m_planes[i].origin();
			HBTK::CartesianLine3D filament(segment.start() + displ, segment.end() + displ);
			// The effect on each of our collocation points.
			for (int n = 0; n < n_planes; n++) {
				wash[n] = wash[n] + 
					HBTK::BiotSavart::unity_vel(collocation_points[n], filament) * vorticity;
			}
		}
	}

	// And now consideration of the streamwise bit of vortex ring (x-dir).
	for (int i = 0; i <= n_planes; i++) {
		double y_pos = (i == 0 ? spanwise_segments[0].start().y() : spanwise_segments[i - 1].end().y());
		HBTK::CartesianLine3D segment = spanwise_segments[i];
		for (int j = 0; j < n_particles_inner; j++) {


			double vorticity = m_vortex_particles[i][j].vorticity;
			for (int n = 0; n < n_planes; n++) {
				wash[n] = wash[n] +
					HBTK::BiotSavart::unity_vel(collocation_points[n], segment) * vorticity;
			}
		}
	}

	return std::vector<double>();
}

std::vector<HBTK::CartesianLine3D> mFlow::LargeAmplitudeLiftingLineTheory::get_segments()
{
	std::vector<HBTK::CartesianLine3D> segments(m_planes.size());
	double start_y = (m_planes[0].origin().y() > 0 ? wing.semispan() : -wing.semispan());
	double end_y = -start_y;
	segments[0].start() = HBTK::CartesianPoint3D({ 0, start_y, 0 });
	segments.back().end() = HBTK::CartesianPoint3D({ 0, end_y, 0 });
	for (int i = 1; i < (int)segments.size(); i++) {
		HBTK::CartesianPoint3D last_orig = m_planes[i-1].origin();
		HBTK::CartesianPoint3D next_orig = m_planes[i].origin();
		HBTK::CartesianPoint3D midpoint = last_orig + (next_orig - last_orig) * 0.5;
		segments[i - 1].end() = midpoint;
		segments[i].start() = midpoint;
	}
	return segments;
}

bool mFlow::LargeAmplitudeLiftingLineTheory::check_matching_vortex_particles()
{
	int n_particles = (int) m_vortex_particles[0].size();
	for (auto & vect : m_vortex_particles) {
		if (n_particles != (int)vect.size()) {
			return false;
		}
	}
	return true;
}
