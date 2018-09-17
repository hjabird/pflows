#include "VortexRingLattice.h"
/*////////////////////////////////////////////////////////////////////////////
VortexRingLattice.cpp

A singular vortex lattice representation.

Copyright 2017-2018 HJA Bird

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

#include <cassert>
#include <cmath>

#include <HBTK/CartesianPlane.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Checks.h>
#include <HBTK/Constants.h>
#include <HBTK/RuntimeProfiler.h>
#include <HBTK/VtkLegacyWriter.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

mFlow::VortexRingLattice::VortexRingLattice(int x_extent, int y_extent)
	: m_geometry(),
	m_ring_strenghts(),
	m_extent_x(x_extent),
	m_extent_y(y_extent)
{
	assert(x_extent >= 0);
	assert(y_extent >= 0);
	m_geometry.extent({ x_extent + 1, y_extent + 1 });
	m_ring_strenghts.extent({ x_extent, y_extent });
}

mFlow::VortexRingLattice::~VortexRingLattice()
{
}

HBTK::CartesianVector3D mFlow::VortexRingLattice::induced_velocity(const HBTK::CartesianPoint3D & measurement_position)
{
	return patch_induced_velocity_inclusive(measurement_position, 0, m_extent_x, 0, m_extent_y);
}

HBTK::CartesianVector3D mFlow::VortexRingLattice::patch_induced_velocity_inclusive(
	const HBTK::CartesianPoint3D & measurement_position,
	int min_x, int max_x, int min_y, int max_y)
{
	HBTK::CartesianVector3D wash({ 0, 0, 0 });
	wash += patch_x_filament_induced_velocity_inclusive(measurement_position,
		min_x, max_x, min_y, max_y);
	wash += patch_y_filament_induced_velocity_inclusive(measurement_position,
		min_x, max_x, min_y, max_y);
	return wash;
}

HBTK::CartesianVector3D mFlow::VortexRingLattice::patch_x_filament_induced_velocity_inclusive(
	const HBTK::CartesianPoint3D & measurement_position,
	int min_x, int max_x, int min_y, int max_y)
{
	assert(min_x >= 0 && min_x <= max_x && max_x <= m_extent_x);
	assert(min_y >= 0 && min_y <= max_y && max_y <= m_extent_y);
	HBTK::CartesianVector3D wash({ 0, 0, 0 });
	for (int ix = min_x; ix < max_x; ix++) {
		for (int iy = min_y; iy <= max_y; iy++) {
			HBTK::CartesianFiniteLine3D filament = edge_x(ix, iy);
			double strength = edge_x_vorticity(ix, iy);
			wash += filament_induced_velocity(filament, strength, measurement_position);
		}
	}
	assert(HBTK::check_finite(wash));
	return wash;
}

HBTK::CartesianVector3D mFlow::VortexRingLattice::patch_y_filament_induced_velocity_inclusive(
	const HBTK::CartesianPoint3D & measurement_position,
	int min_x, int max_x, int min_y, int max_y)
{
	assert(min_x >= 0 && min_x <= max_x && max_x <= m_extent_x);
	assert(min_y >= 0 && min_y <= max_y && max_y <= m_extent_y);
	HBTK::CartesianVector3D wash({ 0, 0, 0 });
	for (int ix = min_x; ix <= max_x; ix++) {
		for (int iy = min_y; iy < max_y; iy++) {
			HBTK::CartesianFiniteLine3D filament = edge_y(ix, iy);
			double strength = edge_y_vorticity(ix, iy);
			wash += filament_induced_velocity(filament, strength, measurement_position);
		}
	}
	assert(HBTK::check_finite(wash));
	return wash;
}

HBTK::CartesianFiniteLine3D mFlow::VortexRingLattice::edge_x(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y <= m_extent_y);
	HBTK::CartesianFiniteLine3D edge;
	edge.start() = m_geometry.value({ index_x, index_y });
	edge.end() = m_geometry.value({ index_x + 1, index_y });
	return edge;
}

double mFlow::VortexRingLattice::edge_x_vorticity(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y <= m_extent_y);
	double vort = 0;
	vort -= (index_y < m_extent_y ? m_ring_strenghts.value({ index_x, index_y }) : 0);
	vort += (index_y > 0 ? m_ring_strenghts.value({ index_x , index_y - 1 }) : 0);
	return vort;
}

HBTK::CartesianFiniteLine3D mFlow::VortexRingLattice::edge_y(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x <= m_extent_x);
	assert(index_y < m_extent_y);
	HBTK::CartesianFiniteLine3D edge;
	edge.start() = m_geometry.value({ index_x, index_y });
	edge.end() = m_geometry.value({ index_x, index_y + 1 });
	return edge;
}

double mFlow::VortexRingLattice::edge_y_vorticity(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x <= m_extent_x);
	assert(index_y < m_extent_y);
	double vort = 0;
	vort -= (index_x < m_extent_x ? m_ring_strenghts.value({ index_x, index_y }) : 0);
	vort += (index_x > 0 ? m_ring_strenghts.value({ index_x - 1, index_y }) : 0);
	return vort;
}

double mFlow::VortexRingLattice::ring_strength(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y < m_extent_y);
	return m_ring_strenghts.value({ index_x, index_y });
}

void mFlow::VortexRingLattice::ring_strength(int index_x, int index_y, double vorticity)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y < m_extent_y);
	m_ring_strenghts.value({ index_x, index_y }) = vorticity;
	return;
}

HBTK::CartesianPoint3D mFlow::VortexRingLattice::vertex(int index_x, int index_y)
{
	return m_geometry.value({ index_x, index_y });
}

void mFlow::VortexRingLattice::vertex(int index_x, int index_y, HBTK::CartesianPoint3D position)
{
	m_geometry.value({ index_x, index_y }) = position;
}

std::array<int, 2> mFlow::VortexRingLattice::extent() const
{
	return std::array<int, 2>({ m_extent_x, m_extent_y });
}

int mFlow::VortexRingLattice::size() const
{
	return m_extent_x * m_extent_y;
}

HBTK::Vtk::VtkUnstructuredDataset mFlow::VortexRingLattice::to_vtk_data()
{
	HBTK::Vtk::VtkUnstructuredDataset data;
	int point_count = 0;
	for (int ix = 0; ix < m_extent_x; ix++) {
		for (int iy = 0; iy <= m_extent_y; iy++) {
			HBTK::CartesianFiniteLine3D edge = edge_x(ix, iy);
			data.mesh.points.push_back(edge.start());
			data.mesh.points.push_back(edge.end());
			data.mesh.cells.push_back({ 3,{ point_count, point_count + 1 } });
			point_count += 2;
			data.scalar_cell_data["Vorticity"].push_back(edge_x_vorticity(ix, iy));
		}
	}
	for (int ix = 0; ix <= m_extent_x; ix++) {
		for (int iy = 0; iy < m_extent_y; iy++) {
			HBTK::CartesianFiniteLine3D edge = edge_y(ix, iy);
			data.mesh.points.push_back(edge.start());
			data.mesh.points.push_back(edge.end());
			data.mesh.cells.push_back({ 3,{ point_count, point_count + 1 } });
			point_count += 2;
			data.scalar_cell_data["Vorticity"].push_back(edge_y_vorticity(ix, iy));
		}
	}
	return data;
}

HBTK::CartesianVector3D mFlow::VortexRingLattice::filament_induced_velocity(
	const HBTK::CartesianFiniteLine3D & filament,
	double strength,
	const HBTK::CartesianPoint3D & measurement_point)
{
	HBTK::CartesianVector3D r1, r2;
	r1 = measurement_point - filament.start();
	r2 = measurement_point - filament.end();
	HBTK::CartesianVector3D r0 = r1 - r2;

	double term_1, term_2, term_21, term_22;
	term_1 = strength / (4 * HBTK::Constants::pi() * (pow(r1.cross(r2).magnitude(), 2)));
	term_21 = r1.dot(r0) / r1.magnitude();
	term_22 = r2.dot(r0) / r2.magnitude();
	term_2 = term_21 - term_22;
	return term_1 * term_2 * r1.cross(r2);
}

