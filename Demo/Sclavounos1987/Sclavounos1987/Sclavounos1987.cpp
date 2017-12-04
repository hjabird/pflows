// Sclavounos1987.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <functional>
#include <iostream>
#include <iomanip>

#include "../../../pFlow/Sclavounos1987.h"
#include "../../../../HBTK/HBTK/GnuPlot.h"
#include "../../../../HBTK/HBTK/Constants.h"


int main()
{
	std::cout << "Starting Sclavounos1987 example.\n\n";

	mFlow::Sclavounos1987 analysis;

	double max_chord = 1;
	analysis.wing.span = 2;
	analysis.wing.m_LE_expr = [&](double x) {
		return -(max_chord / analysis.wing.span)*sqrt(pow(analysis.wing.semispan(), 2) - pow(x,2)); };
	analysis.wing.m_TE_expr = [&](double x) {
		return (max_chord / analysis.wing.span)*sqrt(pow(analysis.wing.semispan(), 2) - pow(x,2)); };
	std::vector<double> X_wing, LE_wing, TE_wing;
	int wing_points = 100;
	X_wing.resize(wing_points); LE_wing.resize(wing_points); TE_wing.resize(wing_points);
	for (int idx = 0; idx < wing_points; idx++) {
		X_wing[idx] = idx * (analysis.wing.span / (wing_points-1)) - analysis.wing.semispan();
		LE_wing[idx] = analysis.wing.leading_edge_X(X_wing[idx]);
		TE_wing[idx] = analysis.wing.trailing_edge_X(X_wing[idx]);
	}
	HBTK::GnuPlot wing_plt;
	wing_plt.hold_on();
	wing_plt.axis_equal_on();
	wing_plt.title("Wing shape projection on XY plane");
	wing_plt.plot(X_wing, LE_wing, "k-");
	wing_plt.plot(X_wing, TE_wing, "k-");

	// Important settings! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	analysis.U = 1;
	analysis.omega = 0;
	analysis.j = 3;
	analysis.number_of_terms = 1;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Elliptic added mass.
	auto added_mass = [](mFlow::Sclavounos1987 analysis)->double {
		return HBTK::Constants::pi() *
			pow(pow(analysis.wing.semispan(), 2) - pow(analysis.wing.semichord(0.), 2), 2)
			/  (4 * analysis.omega);
	};

	
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
	//std::cout << std::setw(10) << "Span pos";
	//std::cout << std::setw(10) << "Chord";
	//std::cout << std::setw(30) << "Circulation\n";

	std::vector<double> X, Y_R, Y_IMG, CL_Span;
	const int n_span_pts = 151;
	X.resize(n_span_pts);
	Y_R.resize(n_span_pts);
	Y_IMG.resize(n_span_pts);
	CL_Span.resize(n_span_pts);
	for (int idx = 0; idx < n_span_pts; idx++) {
		double y = analysis.wing.span * (-0.5 + (double) idx / (n_span_pts-1));
		//std::cout << std::setw(10) << y;
		//std::cout << std::setw(10) << analysis.wing.chord(y);
		//std::cout << std::setw(30) << analysis.get_solution_vorticity(y);
		//std::cout << "\n";
		X[idx] = y;
		Y_R[idx] = analysis.get_solution_vorticity(y).real();
		Y_IMG[idx] = analysis.get_solution_vorticity(y).imag();
		//CL_Span[idx] = analysis.F(y).real() * analysis.wing.semispan();
	}
	std::cout << analysis.m_collocation_points << std::endl;
	HBTK::GnuPlot const_osc_plt, const_osc_F;

	const_osc_plt.replot_off();
	const_osc_plt.plot(X, Y_IMG, ":r");
	const_osc_plt.hold_on();
	const_osc_plt.grid_on();
	const_osc_plt.title("Circulation against span position");
	const_osc_plt.ylabel("Circulation");
	const_osc_plt.xlabel("Span position");
	const_osc_plt.legend({ "Imaginary", "Real" });
	const_osc_plt.plot(X, Y_R, ":b");
	const_osc_plt.replot();

	//const_osc_F.plot(X, CL_Span);

	return 0; // STOP

	const int n_cl_pts = 37;
	const double max_fq_term = 8;
	std::vector<double> abs_cl, phase_res, fq_term;
	abs_cl.resize(n_cl_pts);
	fq_term.resize(n_cl_pts);
	phase_res.resize(n_cl_pts);
#pragma omp parallel for
	for (int idx = 0; idx < n_cl_pts; idx++) {
		mFlow::Sclavounos1987 this_analysis = analysis;
		std::stringstream tmp;
		tmp << "Computing " << idx + 1 << " of " << n_cl_pts << "\n";
		std::cout << tmp.str();
		double tmp_omega = 1e-9 + max_fq_term * (double)idx / (n_cl_pts - 1) * (analysis.U / analysis.wing.semispan());
		this_analysis.omega = tmp_omega;
		this_analysis.compute_solution();
		auto C_L = this_analysis.compute_lift_coeff(added_mass(this_analysis));
		abs_cl[idx] = abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] = tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig7, fig8;
	fig7.title("Figure 5: Modulus of the heave lift-coeffient for elliptical wing of AR = 4");
	fig7.xlabel("wd/U");
	fig7.ylabel("abs(C_L)");
	fig7.plot(fq_term, abs_cl);
	fig8.title("Figure 6: Phase of heave lift coefficeint for elliptical wing of AR = 4");
	fig8.xlabel("wd/U");
	fig8.ylabel("ph(C_L)");
	fig8.plot(fq_term, phase_res);


	system("pause");
}

