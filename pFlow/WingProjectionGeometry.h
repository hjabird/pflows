#pragma once
/*////////////////////////////////////////////////////////////////////////////
WingProjectionGeometry.h

An object to represent a wing projected onto a plane, with geometry 
represented by leading edge and trailing edge functions. 

Copyright 2017 HJA Bird

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

#include <functional>

class WingProjectionGeometry
{
public:
	WingProjectionGeometry();
	~WingProjectionGeometry();

	double leading_edge_X(const double & Y_Global);
	double trailing_edge_X(const double & Y_Global);
	double chord_length(const double & Y_Global);
	double midchord_X(const double & Y_Global);
	
	double cos_angle_between_midchord_and_edge(double Y_global);
	double radius_of_midchord(double Y_global);
	
	double wing_span;
	double standard_chord;

	// Functions c_l(Y) describing the x pos of the leading
	// and trailing edge in terms of y. Assumed valid for 
	// Y = [-B/2, B/2] where B is span.
	std::function<double(const double &)> m_LE_expr;
	std::function<double(const double &)> m_TE_expr;

	// Calculate the std chord based on wing are and span.
	void calculate_standard_chord();
};

