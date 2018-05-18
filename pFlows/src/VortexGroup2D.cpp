#include "VortexGroup2D.h"
/*////////////////////////////////////////////////////////////////////////////
VortexGroup2D.cpp

A class standardise the representation of a set of vortices in 2D from a
single source.

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
#include <cassert>

#include <HBTK/CartesianVector.h>
#include <HBTK/VtkUnstructuredDataset.h>
#include <HBTK/VtkWriter.h>

int mFlow::VortexGroup2D::size() const
{
	return (int) m_vortices.size();
}

const mFlow::VortexGroup2D::Vortex2D & mFlow::VortexGroup2D::operator[](int index) const
{
	return m_vortices[index];
}

mFlow::VortexGroup2D::Vortex2D & mFlow::VortexGroup2D::operator[](int index)
{
	return m_vortices[index];
}

const mFlow::VortexGroup2D::Vortex2D & mFlow::VortexGroup2D::most_recently_added() const
{
	assert(m_last_added_vortex >= 0);	// -1 if invalid.
	return m_vortices[m_last_added_vortex];
}

mFlow::VortexGroup2D::Vortex2D & mFlow::VortexGroup2D::most_recently_added()
{
	assert(m_last_added_vortex >= 0);	// -1 if invalid.
	return m_vortices[m_last_added_vortex];
}

void mFlow::VortexGroup2D::add(const Vortex2D & new_vortex)
{
	m_vortices.emplace_back(new_vortex);
	m_last_added_vortex = (int)m_vortices.size() - 1;
	return;
}

void mFlow::VortexGroup2D::add(const HBTK::CartesianPoint2D & position, const double vorticity)
{
	Vortex2D vort;
	vort.position = position;
	vort.vorticity = vorticity;
	m_vortices.emplace_back(vort);
	m_last_added_vortex = (int)m_vortices.size() - 1;
	return;
}

void mFlow::VortexGroup2D::erase(int index)
{
	m_vortices.erase(m_vortices.begin() + index);
	if (index == m_last_added_vortex) {
		m_last_added_vortex = -1;
	}
	else {
		m_last_added_vortex -= 1;
	}
	return;
}

void mFlow::VortexGroup2D::save_to_vtk(std::ostream & ostream, const HBTK::CartesianPlane & plane) const
{
	HBTK::Vtk::VtkWriter writer;
	HBTK::Vtk::VtkUnstructuredDataset data;
	for (int i = 0; i < size(); i++) {
		auto particle = operator[](i);
		data.mesh.points.push_back(plane(particle.position));
		data.mesh.cells.push_back({ 1, std::vector<int>({ i }) });
		data.scalar_point_data["Vorticity"].push_back(particle.vorticity);
	}
	writer.appended = false;
	writer.open_file(ostream, HBTK::Vtk::VtkWriter::vtk_file_type::UnstructuredGrid);
	writer.write_piece(ostream, data);
	writer.close_file(ostream);
	return;
}

void mFlow::remove_vortices_by_group(std::vector<VortexGroup2D*> groups, const std::vector<std::vector<int>>& indexes_to_remove_by_group)
{
	for (auto & grp_ptr : groups) assert(grp_ptr != NULL);
	assert(groups.size() == indexes_to_remove_by_group.size());
	if (groups.size() == 0) return;

	std::vector<std::vector<int>> sorted_to_remove(indexes_to_remove_by_group);
	for (auto & vect : sorted_to_remove) std::sort(vect.begin(), vect.end());

	int smallest_grp;
	std::vector<int> to_remove;
	{
		int smallest_group_size;
		smallest_grp = 0;
		smallest_group_size = (int)indexes_to_remove_by_group[0].size();
		for (int i = 1; i < (int)indexes_to_remove_by_group.size(); i++) {
			if ((int)indexes_to_remove_by_group[i].size() < smallest_group_size) {
				smallest_group_size = (int)indexes_to_remove_by_group[i].size();
				smallest_grp = i;
			}
		}
	}

	std::vector<int> group_positions;
	for (int & i : group_positions) i = 0;

	int correction = 0;
	for (int i : sorted_to_remove[smallest_grp]) {
		// So sue me for looking searching for a number we took from a group...
		bool remove = true;
		for (auto grp : groups) {
			if (grp->size() <= i - correction) remove = false;
		}
		for (int j = 0; j < (int)groups.size(); j++) {
			while (sorted_to_remove[j][group_positions[j]] < i) {
				group_positions[j]++;
			}
		}
		for (int j = 0; j < (int)groups.size(); j++) {
			if (sorted_to_remove[j][group_positions[j]] != i) {
				remove = false;
			}
		}
		if (remove) {
			for (auto grp : groups) grp->erase(i - correction);
		}
	}
}

void mFlow::remove_vortices_by_group(VortexGroup2D & grp, const std::vector<int>& to_remove)
{
	std::vector<int> remove_sorted(to_remove);
	std::sort(remove_sorted.begin(), remove_sorted.end());
	int correction = 0;
	for (int i : remove_sorted) {
		grp.erase(i - correction);
		correction++;
	}
	return;
}

std::vector<int> mFlow::vorticies_by_critical_distance(VortexGroup2D & group,
	HBTK::CartesianPoint2D reference_point, double critical_distance)
{
	assert(critical_distance >= 0);
	std::vector<int> outside_criterion;
	for (int i = 0; i < group.size(); i++) {
		if ((reference_point - group[i].position).magnitude() > critical_distance) {
			outside_criterion.push_back(i);
		}
	}
	return outside_criterion;
}
