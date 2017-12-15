#include "stdafx.h"
#include "WingProjectionGeometry.h"
/*////////////////////////////////////////////////////////////////////////////
WingProjectionGeometry.cpp

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

#include <cmath>
#include <functional>
#include <cassert>


#include "../../HBTK/HBTK/Integrators.h"
#include "../../HBTK/HBTK/NumericalDifferentiation.h"
#include "../../HBTK/HBTK/Checks.h"

WingProjectionGeometry::WingProjectionGeometry()
{
	invalidate_calculations();
}


WingProjectionGeometry::~WingProjectionGeometry()
{
}


double WingProjectionGeometry::leading_edge_X(const double & Y_coord)
{
	assert(m_LE_expr);
	auto res = m_LE_expr(Y_coord);
	assert(HBTK::check_finite(res));
	return res;
}


double WingProjectionGeometry::trailing_edge_X(const double & Y_coord)
{
	assert(m_TE_expr);
	auto res = m_TE_expr(Y_coord);
	assert(HBTK::check_finite(res));
	return res;
}

double WingProjectionGeometry::chord(const double & Y_Global)
{
	assert(m_LE_expr);
	assert(m_TE_expr);
	double TE_x = trailing_edge_X(Y_Global);
	double LE_x = leading_edge_X(Y_Global);
	assert(TE_x >= LE_x);
	return TE_x - LE_x;
}

double WingProjectionGeometry::semichord(const double & Y_Global)
{
	return chord(Y_Global) / 2;
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
	double dydx = HBTK::central_difference_O1A2(fn, Y_global);
	assert(HBTK::check_finite(dydx));
	return cos(atan(dydx));
}

double WingProjectionGeometry::radius_of_midchord(double Y_global)
{
	std::function<double(double)> fn = [&](double y) {return midchord_X(y); };
	double dx_term = HBTK::central_difference_O1A2(fn, Y_global);
	double ddx_term = HBTK::central_difference_O2A2(fn, Y_global);

	return abs(pow(sqrt(1 + dx_term), 3) / ddx_term);
}

const double WingProjectionGeometry::semispan()
{
	return span / 2;
}

double WingProjectionGeometry::standard_chord()
{
	if (!m_valid_standard_chord) { 
		calculate_standard_chord();
		m_valid_standard_chord = true;
	}
	return m_standard_chord;
}

double WingProjectionGeometry::area()
{
	if (!m_valid_area) {
		calculate_wing_area();
		m_valid_area = true;
	}
	return m_area;
}

double WingProjectionGeometry::aspect_ratio()
{
	if (!m_valid_aspect_ratio)
	{
		calculate_aspect_ratio();
	}
	return m_aspect_ratio;
}

void WingProjectionGeometry::add_to_plot(HBTK::GnuPlot &plot) {
	plot.hold_on();
	plot.replot_off();
	plot.plot([&](double x)->double {return m_LE_expr(x); }, -semispan(), semispan(), "k-");
	plot.plot([&](double x)->double {return m_TE_expr(x); }, -semispan(), semispan(), "k-");
	plot.replot_on();
	plot.axis_equal_on();
	plot.hold_off();
	plot.replot();
	return;
}

void WingProjectionGeometry::invalidate_calculations()
{
	m_valid_area = false;
	m_valid_standard_chord = false;
	m_valid_aspect_ratio = false;
}

void WingProjectionGeometry::calculate_standard_chord()
{
	m_standard_chord = area() / span;
	return;
}

void WingProjectionGeometry::calculate_wing_area()
{
	assert(m_LE_expr);
	assert(m_TE_expr);
	auto my_function = [&](const double y_loc) {
		auto c = abs(trailing_edge_X(y_loc) - leading_edge_X(y_loc));
		return c;
	};
	m_area = HBTK::adaptive_simpsons_integrate(my_function, 1e-9, -semispan(), semispan());
	assert(HBTK::check_finite(m_area));
	return;
}

void WingProjectionGeometry::calculate_aspect_ratio()
{
	m_aspect_ratio = span*span / area();
	m_valid_aspect_ratio = true;
	return;
}
