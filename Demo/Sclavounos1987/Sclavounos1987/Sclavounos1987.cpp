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

	// Elliptic
	analysis.wing.span = 4.;
	auto aspect_ratio = 4;
	auto max_chord = 4 * analysis.wing.span / (aspect_ratio * HBTK::Constants::pi());
	analysis.wing.m_LE_expr = [&](double x) {
		return -(max_chord / analysis.wing.span)*sqrt(pow(analysis.wing.semispan(), 2) - pow(x, 2)); };
	analysis.wing.m_TE_expr = [&](double x) {
		return (max_chord / analysis.wing.span)*sqrt(pow(analysis.wing.semispan(), 2) - pow(x, 2)); };
	/*
	analysis.wing.span = 4;
	analysis.wing.m_LE_expr = [&](double x) { return -0.5; };
	analysis.wing.m_TE_expr = [&](double x) { return 0.5; };
	*/
	std::vector<double> X_wing, LE_wing, TE_wing;
	int wing_points = 100;
	X_wing.resize(wing_points); LE_wing.resize(wing_points); TE_wing.resize(wing_points);
	for (int idx = 0; idx < wing_points; idx++) {
		X_wing[idx] = idx * (analysis.wing.span / (wing_points - 1)) - analysis.wing.semispan();
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
	analysis.omega = 8;
	analysis.j = 3;
	analysis.number_of_terms = 8;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Elliptic added mass.
	auto added_mass = [](mFlow::Sclavounos1987 analysis)->double {
		// Reference: http://brennen.caltech.edu/fluidbook/basicfluiddynamics/unsteadyflows/addedmass/valuesoftheaddedmass.pdf
		// (Unable to find any closed form solution (even in Hydrodynamics(Lamb) 2nd Ed.)
		// Integration of the flat plate term in eq4.4 results in the equivalent of AR=1, without the correction
		// for the elliptic nature of the plate.
		auto a = analysis.wing.semispan();
		auto b = analysis.wing.semichord(0.);
		if (a < b) { std::swap(a, b); }
		double added_mass = a * b * b * HBTK::Constants::pi() * 4. / 3.;
		double ratio = a / b;
		std::vector<double> known_ratios = { 1., 1.5, 2., 3., 4., 6., 8.19, 10.34, 14.30, 10000. }; // To inf really.
		std::vector<double> known_coeffs = { 0.637, 0.748, 0.826, 0.900, 0.933, 0.964, 0.978, 0.985, 0.991, 1.000 };

		double coeff = 0;
		int lower_known = -1;
		// Linear interpolation
		for (int idx = 0; idx < (int)known_ratios.size(); idx++) {
			if (ratio >= known_ratios[idx]) {
				lower_known = idx;
			}
		}
		if (lower_known == 0) { coeff = known_ratios[0]; } // == 1.

		double fraction = (ratio - known_ratios[lower_known]) / (known_ratios[lower_known + 1] - known_ratios[lower_known]);
		coeff = (known_coeffs[lower_known + 1] - known_coeffs[lower_known]) * fraction + known_coeffs[lower_known];
		added_mass *= 2 * coeff;
		return added_mass;
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

	std::vector<double> X, Y_R, Y_IMG;
	const int n_span_pts = 151;
	X.resize(n_span_pts);
	Y_R.resize(n_span_pts);
	Y_IMG.resize(n_span_pts);
	for (int idx = 0; idx < n_span_pts; idx++) {
		double y = analysis.wing.span * (-0.5 + (double)idx / (n_span_pts - 1));
		X[idx] = y;
		Y_R[idx] = analysis.solution_vorticity(y).real();
		Y_IMG[idx] = analysis.solution_vorticity(y).imag();
	}

	HBTK::GnuPlot const_osc_plt, const_osc_F;

	const_osc_plt.replot_off();
	const_osc_plt.plot(X, Y_IMG, ":r");
	const_osc_plt.hold_on();
	const_osc_plt.grid_on();
	const_osc_plt.title("Circulation against span position (elliptic wing AR=" + std::to_string(analysis.wing.aspect_ratio()) +
		", omega=" + std::to_string(analysis.omega) + ", U=" + std::to_string(analysis.U));
	const_osc_plt.ylabel("Circulation");
	const_osc_plt.xlabel("Span position");
	const_osc_plt.legend({ "Imaginary", "Real" });
	const_osc_plt.plot(X, Y_R, ":b");
	const_osc_plt.replot();

	const int n_cl_pts = 150;
	const double max_fq_term = 8;
	std::vector<double> abs_cl, phase_res, fq_term, cl_i, cl_r;
	abs_cl.resize(n_cl_pts);
	fq_term.resize(n_cl_pts);
	phase_res.resize(n_cl_pts);
	cl_i.resize(n_cl_pts);
	cl_r.resize(n_cl_pts);

#pragma omp parallel for
	for (int idx = 0; idx < n_cl_pts; idx++) {
		mFlow::Sclavounos1987 this_analysis = analysis;
		std::stringstream tmp;
		tmp << "Computing " << idx + 1 << " of " << n_cl_pts << "\n";
		std::cout << tmp.str();
		double tmp_omega = 1e-14 + max_fq_term * (double)idx / (n_cl_pts - 1) * (analysis.U / analysis.wing.semispan());
		this_analysis.omega = tmp_omega;
		this_analysis.compute_solution();
		auto C_L = this_analysis.compute_lift_coeff(added_mass(this_analysis));
		cl_r[idx] = C_L.real();
		cl_i[idx] = C_L.imag();
		abs_cl[idx] = abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] = tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig7, fig8;
	fig7.title("Figure 5: Modulus of the heave lift-coeffient for elliptical wing of AR = 4");
	fig7.xlabel("wd/U");
	fig7.ylabel("abs(C_L)");
	fig7.hold_on();
	fig7.plot(fq_term, abs_cl, "k-");
	fig7.plot(fq_term, cl_r, "r-");
	fig7.plot(fq_term, cl_i, "b-");
	fig7.legend({ "abs(Cl)", "Real(Cl)", "Imag(Cl)" });
	fig7.replot();

	fig8.title("Figure 6: Phase of heave lift coefficeint for elliptical wing of AR = 4");
	fig8.xlabel("wd/U");
	fig8.ylabel("ph(C_L)");
	fig8.plot(fq_term, phase_res, "k-");


	return 0;
}
