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

#include <vector>

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

	private:
		std::vector<Vortex2D> m_vortices;	// 
		int m_last_added_vortex;		// = -1 for invalid.
	};
}