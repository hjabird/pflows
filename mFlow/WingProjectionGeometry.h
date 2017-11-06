#pragma once
/*////////////////////////////////////////////////////////////////////////////
WingProjectionGeometry.h

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

