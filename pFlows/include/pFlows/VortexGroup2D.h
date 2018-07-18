#pragma once
/*////////////////////////////////////////////////////////////////////////////
VortexGroup2D.h

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

#include <ostream>
#include <vector>

#include <HBTK/CartesianPlane.h>
#include <HBTK/CartesianPoint.h>

namespace mFlow {
	class VortexGroup2D {
	public:
		struct Vortex2D {
			HBTK::CartesianPoint2D position;
			double vorticity;
		};

		// Number of vortices that make up the vortex group.
		int size() const;

		// Total vorticity of the group
		double vorticity_sum() const;

		// Access a vortex (index in [0, number_of_vortices) ).
		const Vortex2D & operator[](int index) const;
		Vortex2D & operator[](int index);
		const Vortex2D & most_recently_added() const;
		Vortex2D & most_recently_added();

		// Add a vortex
		void add(const Vortex2D & new_vortex);
		void add(const HBTK::CartesianPoint2D & position, const double vorticity);

		// Remove a vortex by index
		void erase(int index);

		// Save to .vtu file
		void save_to_vtk(std::ostream & ostream, const HBTK::CartesianPlane & plane) const;

	private:
		std::vector<Vortex2D> m_vortices;	// 
		int m_last_added_vortex;		// = -1 for invalid.
	};

	// Remove vortices with corresponding indexes from multiple groups
	// if the criterion is considered in AND fashion.
	void remove_vortices_by_group(std::vector<VortexGroup2D*> groups,
		const std::vector<std::vector<int>> & indexes_to_remove_by_group);

	// Remove vortices with corresponding indexes from group
	void remove_vortices_by_group(VortexGroup2D & grp,
		const std::vector<int> & to_remove);

	// Get the indices of vortices outside some critical distance.
	std::vector<int> vorticies_by_critical_distance(VortexGroup2D & group, 
		HBTK::CartesianPoint2D reference_point, double critical_distance);
}