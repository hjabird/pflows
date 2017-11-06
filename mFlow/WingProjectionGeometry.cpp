#include "stdafx.h"
#include "WingProjectionGeometry.h"
/*////////////////////////////////////////////////////////////////////////////
WingProjectionGeometry.cpp

An object to represent a wing projected onto a plane, with geometry
represented by leading edge and trailing edge functions.

Copyright 2017 HJA Bird

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <functional>
#include <cassert>


#include "../../HBTK/HBTK/Integrators.h"
#include "../../HBTK/HBTK/NumericalDifferentiation.h"

WingProjectionGeometry::WingProjectionGeometry()
{
	// THe values are public accessible, so we can't check they've been initialised.
	// Hence, make them invalid.
	standard_chord = NAN;
	wing_span = NAN;
}


WingProjectionGeometry::~WingProjectionGeometry()
{
}


double WingProjectionGeometry::leading_edge_X(const double & Y_coord)
{
	assert(m_LE_expr);
	return m_LE_expr(Y_coord);
}


double WingProjectionGeometry::trailing_edge_X(const double & Y_coord)
{

	assert(m_TE_expr);
	return m_TE_expr(Y_coord);
}

double WingProjectionGeometry::chord_length(const double & Y_Global)
{
	assert(m_LE_expr);
	assert(m_TE_expr);
	return trailing_edge_X(Y_Global) - leading_edge_X(Y_Global);
}

double WingProjectionGeometry::midchord_X(const double & Y_Global)
{
	assert(m_LE_expr);
	assert(m_TE_expr);
	return (leading_edge_X(Y_Global) + trailing_edge_X(Y_Global)) / 2.0;
}

double WingProjectionGeometry::cos_angle_between_midchord_and_edge(double Y_global)
{
	std::function<double(double)> fn = [&](double y){return midchord_X(y); };
	double dydx = Diff::central_difference_O1A2(fn, Y_global);
	return cos(atan(dydx));
}

double WingProjectionGeometry::radius_of_midchord(double Y_global)
{
	std::function<double(double)> fn = [&](double y) {return midchord_X(y); };
	double dx_term = Diff::central_difference_O1A2(fn, Y_global);
	double ddx_term = Diff::central_difference_O2A2(fn, Y_global);

	return abs(pow(sqrt(1 + dx_term), 3) / ddx_term);
}


void WingProjectionGeometry::calculate_standard_chord()
{
	assert(m_LE_expr);
	assert(m_TE_expr);
	auto my_function = [&](const double y_loc) {
		auto c = leading_edge_X(y_loc) - trailing_edge_X(y_loc);
		return c / 2;
	};
	standard_chord = Quad::adaptive_simpsons_integrate(my_function, 1e-10, -1.0, 1.0);
	return;
}
