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

#include <cassert>

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
