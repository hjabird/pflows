// Sclavounos1987.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <functional>
#include <iostream>
#include <iomanip>

#include "../../../pFlow/Sclavounos1987.h"


int main()
{
	std::cout << "Starting Sclavounos1987 example.\n\n";

	mFlow::Sclavounos1987 analysis;

	analysis.wing.m_LE_expr = [](double x) {return 0.0; };
	analysis.wing.m_TE_expr = [](double x) {return 1.0; };
	analysis.wing.wing_span = 11.;

	analysis.U = 2;
	analysis.omega = 1;
	analysis.j = 3;
	analysis.number_of_terms = 8;
	
	analysis.compute_solution();

	std::cout << "INPUT PARAMATERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	auto pl1 = [&](auto str, auto val)->void
	{
		std::cout << "\t" << std::setw(10) << str << std::setw(15) << val << "\n"; 
	};
	pl1("U", analysis.U);
	pl1("Pert freq.", analysis.omega);
	pl1("j", analysis.j);
	pl1("Num terms", analysis.number_of_terms);

	std::cout << "\n\n";

	std::cout << "RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	std::cout << std::setw(10) << "Span pos";
	std::cout << std::setw(10) << "Chord";
	std::cout << std::setw(30) << "Circulation\n";

	const int n_span_pts = 21;
	for (int idx = 0; idx < n_span_pts; idx++) {
		double y = analysis.wing.wing_span * (-0.5 + (double) idx / (n_span_pts-1));
		std::cout << std::setw(10) << y;
		std::cout << std::setw(10) << analysis.wing.chord_length(y); // This may be wrong?
		std::cout << std::setw(30) << analysis.get_solution_vorticity(y);
		std::cout << "\n";
	}



	system("pause");
}

