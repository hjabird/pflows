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
		double chord = wing.chord(plane_origin.y);
		m_chord_lengths[i] = chord;
	}
}

std::vector<std::vector<HBTK::CartesianRectilinearPanel>> mFlow::LargeAmplitudeLiftingLineTheory::panelise_wake()
{
	std::vector<std::vector<HBTK::CartesianRectilinearPanel>> panels;
	panels.resize(m_planes.size());
	for (auto & vect : panels) vect.resize(m_vortex_particles[0].size() - 1);
	std::vector<double> y_values(panels.size() + 1);
	auto segments = get_segments();
	for (int i = 0; i < panels.size(); i++) {
		y_values[i] = segments[i].start.y;
	}
	y_values.back() = segments.back().end.y;

	// We want to make our mesh as smooth as possible. This requires some work.
	// Thinking of the mesh, so long as the vortex filaments go through all our 
	// points, we're good. Except that actually it is unconstrained and could
	// "zigzag". Here, we try and compute the least squares of the derivatives
	// of the mesh positions.

	// Calculates the i=0 displacement: (Notes #3/pg34)
	auto calc_xm0 = [&](std::vector<double> x_values)->double {
		std::vector<double> y_i(x_values.size());
		for (int i = 0; i < y_i.size(); i++) {
			y_i[i] = y_values[i + 1] - y_values[i];
		}

		double denominator = 0;
		for (double y : y_i) denominator += 2 * 1 / (y*y);

		double numerator = 0;
		for (int i = 0; i < y_i.size(); i++) {
			double internal_sum = 0;
			for (int j = 0; j < i; j++) {
				internal_sum += 2 * x_values[j] * pow(-1, -j - 1);
			}
			numerator += 2 * (x_values[i] - 2 * internal_sum) / (y_i[i] * y_i[i]);
		}
		return numerator / denominator;
	};

	// Panels all start on the lifting line:
	panels[0][0].corners[0] = segments[0].start;
	for (int i = 1; i < panels.size(); i++) {
		panels[i][0].corners[0] = segments[i].start;
		panels[i - 1][0].corners[1] = segments[i].start;
	}
	panels.back()[0].corners[1];

	// Now interate in streamwise direction.
	for (int i = 0; i < m_vortex_particles[0].size(); i++) {
		// Get all the x_i, y_i values (again see notes #3pg34)



		for (int j = 1; j < panels.size(); j++) {
			
		}
	}



	return panels;
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
		// Wakewise vortex ring consideration
		for (int j = 0; j < n_particles_inner; j++) {
			double vorticity = m_vortex_particles[i][j].vorticity;
			// The effect on each of our collocation points.
			for (int n = 0; n < n_planes; n++) {
				wash[n] = wash[n] + 
					HBTK::BiotSavart::unity_vel(collocation_points[n], filament) * vorticity;
			}
		}
	}

	// And now consideration of the streamwise bit of vortex ring (x-dir).
	for (int i = 1; i < n_planes; i++) {
		HBTK::CartesianLine3D segment = spanwise_segments[i];
		// The strength of the adjacent vortex rings must computed.
		double vorticity_sum_ym = 0;
		double vorticity_sum_yp = 0;
		for (int j = 0; j < n_particles_inner; j++) {
			vorticity_sum_ym += m_vortex_particles[i - 1][j].vorticity;
			vorticity_sum_yp += m_vortex_particles[i][j].vorticity;
			double vorticity = vorticity_sum_yp - vorticity_sum_ym;
			
			for (int n = 0; n < n_planes; n++) {
				wash[n] = wash[n] +
					HBTK::BiotSavart::unity_vel(collocation_points[n], segment) * vorticity;
			}
		}
	}

	return std::vector<double>();
}

std::vector<HBTK::CartesianFiniteLine3D> mFlow::LargeAmplitudeLiftingLineTheory::get_segments()
{
	std::vector<HBTK::CartesianFiniteLine3D> segments(m_planes.size());
	double start_y = (m_planes[0].origin().y > 0 ? wing.semispan() : -wing.semispan());
	double end_y = -start_y;
	segments[0].start = HBTK::CartesianPoint3D({ 0, start_y, 0 });
	segments.back().end = HBTK::CartesianPoint3D({ 0, end_y, 0 });
	for (int i = 1; i < (int)segments.size(); i++) {
		HBTK::CartesianPoint3D last_orig = m_planes[i-1].origin();
		HBTK::CartesianPoint3D next_orig = m_planes[i].origin();
		HBTK::CartesianPoint3D midpoint = last_orig + (next_orig - last_orig) * 0.5;
		segments[i - 1].end = midpoint;
		segments[i].start = midpoint;
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
