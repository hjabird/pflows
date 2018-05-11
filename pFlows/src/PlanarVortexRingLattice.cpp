#include "PlanarVortexRingLattice.h"


#include <cassert>

#include <HBTK/CartesianPlane.h>
#include <HBTK/CartesianVector.h>
#include <HBTK/Checks.h>
#include <HBTK/Constants.h>
#include <HBTK/VtkLegacyWriter.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

mFlow::PlanarVortexRingLattice::PlanarVortexRingLattice(int x_extent, int y_extent)
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

mFlow::PlanarVortexRingLattice::~PlanarVortexRingLattice()
{
}

double mFlow::PlanarVortexRingLattice::downwash(const HBTK::CartesianPoint2D & measurement_position)
{
	return patch_downwash_inclusive(measurement_position, 0, m_extent_x, 0, m_extent_y);
}

double mFlow::PlanarVortexRingLattice::patch_downwash_inclusive(
	const HBTK::CartesianPoint2D & measurement_position, 
	int min_x, int max_x, int min_y, int max_y)
{
	double wash = 0;
	wash += patch_x_filament_downwash_inclusive(measurement_position,
		min_x, max_x, min_y, max_y);
	wash += patch_y_filament_downwash_inclusive(measurement_position,
		min_x, max_x, min_y, max_y);
	return wash;
}

double mFlow::PlanarVortexRingLattice::patch_x_filament_downwash_inclusive(
	const HBTK::CartesianPoint2D & measurement_position, 
	int min_x, int max_x, int min_y, int max_y)
{
	assert(min_x >= 0 && min_x <= max_x && max_x <= m_extent_x);
	assert(min_y >= 0 && min_y <= max_y && max_y <= m_extent_y);
	double wash = 0;
	for (int ix = min_x; ix < max_x; ix++) {
		for (int iy = min_y; iy <= max_y; iy++) {
			HBTK::CartesianFiniteLine2D filament = edge_x(ix, iy);
			double strength = edge_x_vorticity(ix, iy);
			wash += filament_downwash(filament, strength, measurement_position);
		}
	}
	assert(HBTK::check_finite(wash));
	return wash;
}

double mFlow::PlanarVortexRingLattice::patch_y_filament_downwash_inclusive(
	const HBTK::CartesianPoint2D & measurement_position, 
	int min_x, int max_x, int min_y, int max_y)
{
	assert(min_x >= 0 && min_x <= max_x && max_x <= m_extent_x);
	assert(min_y >= 0 && min_y <= max_y && max_y <= m_extent_y);
	double wash = 0;
	for (int ix = min_x; ix <= max_x; ix++) {
		for (int iy = min_y; iy < max_y; iy++) {
			HBTK::CartesianFiniteLine2D filament = edge_y(ix, iy);
			double strength = edge_y_vorticity(ix, iy);
			wash += filament_downwash(filament, strength, measurement_position);
		}
	}
	assert(HBTK::check_finite(wash));
	return wash;
}

HBTK::CartesianFiniteLine2D mFlow::PlanarVortexRingLattice::edge_x(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y <= m_extent_y);
	HBTK::CartesianFiniteLine2D edge;
	edge.start() = m_geometry.value({ index_x, index_y });
	edge.end() = m_geometry.value({ index_x + 1, index_y });
	return edge;
}

double mFlow::PlanarVortexRingLattice::edge_x_vorticity(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y <= m_extent_y);
	double vort = 0;
	vort -= (index_y < m_extent_y ? m_ring_strenghts.value({ index_x, index_y })  : 0);
	vort += (index_y > 0 ? m_ring_strenghts.value({ index_x , index_y - 1}) : 0);
	return vort;
}

HBTK::CartesianFiniteLine2D mFlow::PlanarVortexRingLattice::edge_y(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x <= m_extent_x);
	assert(index_y < m_extent_y);
	HBTK::CartesianFiniteLine2D edge;
	edge.start() = m_geometry.value({ index_x, index_y });
	edge.end() = m_geometry.value({ index_x, index_y + 1 });
	return edge;
}

double mFlow::PlanarVortexRingLattice::edge_y_vorticity(int index_x, int index_y)
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

double mFlow::PlanarVortexRingLattice::ring_strength(int index_x, int index_y)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y < m_extent_y);
	return m_ring_strenghts.value({ index_x, index_y });
}

void mFlow::PlanarVortexRingLattice::ring_strength(int index_x, int index_y, double vorticity)
{
	assert(index_x >= 0);
	assert(index_y >= 0);
	assert(index_x < m_extent_x);
	assert(index_y < m_extent_y);
	m_ring_strenghts.value({ index_x, index_y }) = vorticity;
	return;
}

HBTK::CartesianPoint2D mFlow::PlanarVortexRingLattice::vertex(int index_x, int index_y)
{
	return m_geometry.value({ index_x, index_y });
}

void mFlow::PlanarVortexRingLattice::vertex(int index_x, int index_y, HBTK::CartesianPoint2D position)
{
	m_geometry.value({ index_x, index_y }) = position;
}

std::array<int, 2> mFlow::PlanarVortexRingLattice::extent() const
{
	return std::array<int, 2>({ m_extent_x, m_extent_y });
}

int mFlow::PlanarVortexRingLattice::size() const
{
	return m_extent_x * m_extent_y;
}

void mFlow::PlanarVortexRingLattice::save_to_vtk(std::ostream & out_stream)
{
	HBTK::Vtk::VtkWriter writer;
	HBTK::Vtk::VtkUnstructuredDataset data;
	HBTK::CartesianPlane plane(
		HBTK::CartesianPoint3D({ 0, 0, 0 }),
		HBTK::CartesianPoint3D({ 1, 0, 0 }),
		HBTK::CartesianPoint3D({ 0, 1, 0 })
	);
	int point_count = 0;
	for (int ix = 0; ix < m_extent_x; ix++) {
		for (int iy = 0; iy <= m_extent_y; iy++) {
			HBTK::CartesianFiniteLine2D edge = edge_x(ix, iy);
			data.mesh.points.push_back(plane(edge.start()));
			data.mesh.points.push_back(plane(edge.end()));
			data.mesh.cells.push_back({ 3, {point_count, point_count + 1} });
			point_count += 2;
			data.scalar_cell_data["Vorticity"].push_back(edge_x_vorticity(ix, iy));
		}
	}
	for (int ix = 0; ix <= m_extent_x; ix++) {
		for (int iy = 0; iy < m_extent_y; iy++) {
			HBTK::CartesianFiniteLine2D edge = edge_y(ix, iy);
			data.mesh.points.push_back(plane(edge.start()));
			data.mesh.points.push_back(plane(edge.end()));
			data.mesh.cells.push_back({ 3,{ point_count, point_count + 1 } });
			point_count += 2;
			data.scalar_cell_data["Vorticity"].push_back(edge_y_vorticity(ix, iy));
		}
	}
	writer.appended = false;
	writer.open_file(out_stream, HBTK::Vtk::VtkWriter::UnstructuredGrid);
	writer.write_piece(out_stream, data);
	writer.close_file(out_stream);
	return;
}

double mFlow::PlanarVortexRingLattice::filament_downwash(
	const HBTK::CartesianFiniteLine2D & filament, 
	double strength, 
	const HBTK::CartesianPoint2D & measurement_point)
{
	HBTK::CartesianVector2D r1, r2;
	r1 = measurement_point - filament.start();
	r2 = measurement_point - filament.end() ;
	HBTK::CartesianVector2D r0 = r1 - r2;

	double term_1, term_2, term_21, term_22;
	term_1 = strength / (4 * HBTK::Constants::pi() * (r1.x() * r2.y() - r1.y() * r2.x()));
	term_21 = r1.dot(r0) / r1.magnitude();
	term_22 = r2.dot(r0) / r2.magnitude();
	term_2 = term_21 - term_22;
	return term_1 * term_2;
}

