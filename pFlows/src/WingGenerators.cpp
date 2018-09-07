#include "WingGenerators.h"
/*////////////////////////////////////////////////////////////////////////////
WingGenerators.cpp

Methods to generate wing shapes.

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

#include <HBTK/Checks.h>
#include <HBTK/Constants.h>

void mFlow::WingGenerators::elliptic(WingProjectionGeometry & wing_obj, double span, double aspect_ratio)
{
	span = std::abs(span);
	aspect_ratio = std::abs(aspect_ratio);
	auto max_chord = 4 * span / (aspect_ratio * HBTK::Constants::pi());
	auto leading_edge = [=](double x_position)->double {
		double p =  -(max_chord / span)*sqrt(pow(span / 2, 2) - pow(x_position, 2));
		if (!HBTK::check_finite(p)) p = 0.0;
		return p;
	};
	auto trailing_edge = [=](double x_position)->double {
		double p = (max_chord / span)*sqrt(pow(span / 2, 2) - pow(x_position, 2));		
		if (!HBTK::check_finite(p)) p = 0.0;
		return p;
	};

	wing_obj.span = span;
	wing_obj.m_TE_expr = trailing_edge;
	wing_obj.m_LE_expr = leading_edge;

	return;
}

void mFlow::WingGenerators::rectangular(WingProjectionGeometry & wing_obj, double span, double aspect_ratio)
{
	span = std::abs(span);
	aspect_ratio = std::abs(aspect_ratio);
	auto max_chord = span / aspect_ratio;
	auto leading_edge = [=](double x_position)->double {
		return -max_chord / 2;
	};
	auto trailing_edge = [=](double x_position)->double {
		return max_chord / 2;
	};

	wing_obj.span = span;
	wing_obj.m_TE_expr = trailing_edge;
	wing_obj.m_LE_expr = leading_edge;

	return;
}

