// Sclavounos1987.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"

#include <functional>
#include <iostream>

#include "../../../pFlow/Sclavounos1987.h"
#include "../../../pFlow/WingGenerators.h"
#include "../../../../HBTK/HBTK/GnuPlot.h"
#include "../../../../HBTK/HBTK/Constants.h"

void run_elliptic_experiment() {

	mFlow::Sclavounos1987 analysis;

	mFlow::WingGenerators::elliptic(analysis.wing, 4., 4.);
	analysis.U = 1;
	analysis.omega = 0.001;
	analysis.j = 3;
	analysis.number_of_terms = 8;

	std::cout << "\n\n";
	std::cout << "Computing results for comparison to Sclavounos1987 figures 5 & 6.\n";
	std::cout << "Using elliptic wing of AR = " << analysis.wing.aspect_ratio() << 
		" and span = " << analysis.wing.span << ".\n";
	std::cout << "Using " << analysis.number_of_terms << " fourier terms.\n";
	std::cout << "Using j = " << analysis.j << ".\n";
	std::cout << "Using U = " << analysis.U << ".\n\n";

	const int number_of_points = 150;
	const double max_fq_term = 8;
	std::vector<double> cl_abs, phase_res, fq_term, cl_imag, cl_real;
	cl_abs.resize(number_of_points);
	fq_term.resize(number_of_points);
	phase_res.resize(number_of_points);
	cl_imag.resize(number_of_points);
	cl_real.resize(number_of_points);

#pragma omp parallel for
	for (int idx = 0; idx < number_of_points; idx++) {
		mFlow::Sclavounos1987 this_analysis = analysis;
		double tmp_omega = 1e-14 + max_fq_term * (double)idx / (number_of_points - 1) * (analysis.U / analysis.wing.semispan());
		this_analysis.omega = tmp_omega;
		this_analysis.compute_solution();
		std::complex<double> C_L = this_analysis.compute_lift_coeff_j3(this_analysis.elliptic_added_mass_coefficient());
		cl_real[idx] =	C_L.real();
		cl_imag[idx] =	C_L.imag();
		cl_abs[idx] =	abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] =	tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig5, fig6;
	fig5.title("Figure 5: Modulus of the heave lift-coeffient for elliptical wing of AR = 4");
	fig5.xlabel("wd/U");
	fig5.ylabel("abs(C_L)");
	fig5.yrange(0.0, 8.0);
	fig5.hold_on();
	fig5.plot(fq_term, cl_abs, "k-");
	fig5.plot(fq_term, cl_real, "r-");
	fig5.plot(fq_term, cl_imag, "b-");
	fig5.legend({ "abs(Cl)", "Real(Cl)", "Imag(Cl)" });
	fig5.replot();

	fig6.title("Figure 6: Phase of heave lift coefficeint for elliptical wing of AR = 4");
	fig6.xlabel("wd/U");
	fig6.ylabel("ph(C_L) (deg)");
	fig6.yrange(-270, 0);
	fig6.plot(fq_term, phase_res, "k-");

	system("pause");
	return;
} 

void run_rectangular_experiment() {

	mFlow::Sclavounos1987 analysis;
	mFlow::WingGenerators::rectangular(analysis.wing, 4, 4);
	analysis.U = 1;
	analysis.omega = 0.001;
	analysis.j = 3;
	analysis.number_of_terms = 8;

	std::cout << "\n\n";
	std::cout << "Computing results for comparison to Sclavounos1987 figures 7 & 8.\n";
	std::cout << "Using rectangular wing of AR = " << analysis.wing.aspect_ratio() <<
		" and span = " << analysis.wing.span << ".\n";
	std::cout << "Using " << analysis.number_of_terms << " fourier terms.\n";
	std::cout << "Using j = " << analysis.j << ".\n";
	std::cout << "Using U = " << analysis.U << ".\n\n";

	const int number_of_points = 150;
	const double max_fq_term = 8;
	std::vector<double> cl_abs, phase_res, fq_term, cl_imag, cl_real;
	cl_abs.resize(number_of_points);
	fq_term.resize(number_of_points);
	phase_res.resize(number_of_points);
	cl_imag.resize(number_of_points);
	cl_real.resize(number_of_points);

#pragma omp parallel for
	for (int idx = 0; idx < number_of_points; idx++) {
		mFlow::Sclavounos1987 this_analysis = analysis;
		double tmp_omega = 1e-14 + max_fq_term * (double)idx / (number_of_points - 1) * (analysis.U / analysis.wing.semispan());
		this_analysis.omega = tmp_omega;
		this_analysis.compute_solution();
		std::complex<double> C_L = 
			this_analysis.compute_lift_coeff_j3(this_analysis.rectangular_added_mass_coefficient());
		cl_real[idx] = C_L.real();
		cl_imag[idx] = C_L.imag();
		cl_abs[idx] = abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] = tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig7, fig8;
	fig7.title("Figure 7: Modulus of the heave lift-coeffient for rectangular wing of AR = 4");
	fig7.xlabel("wd/U");
	fig7.ylabel("abs(C_L)");
	fig7.yrange(0, 8);
	fig7.hold_on();
	fig7.plot(fq_term, cl_abs, "k-");
	fig7.plot(fq_term, cl_real, "r-");
	fig7.plot(fq_term, cl_imag, "b-");
	fig7.legend({ "abs(Cl)", "Real(Cl)", "Imag(Cl)" });
	fig7.replot();

	fig8.title("Figure 8: Phase of heave lift coefficeint for rectangular wing of AR = 4");
	fig8.xlabel("wd/U");
	fig8.ylabel("ph(C_L)");
	fig8.yrange(-270, 0);
	fig8.plot(fq_term, phase_res, "k-");

	system("pause");
	return;
}


int main()
{
	std::cout << "Starting Sclavounos1987 example.\n\n";
	std::cout << "This code is based on An unsteady lifting-line theory, P.D. Sclavounos, 1987,\n";
	std::cout << "Journal of Engineering Mathmatics. https://doi.org/10.1007/BF00127464\n";
	std::cout << "Code copyright HJA Bird 2017. Available under GNU GPL3 lisence from \n";
	std::cout << "https://github.com/hjabird/pflows. For this to work it is required that\n";
	std::cout << "GnuPlot is available and added to path.\n\n";

	run_elliptic_experiment();
	run_rectangular_experiment();

	return 0;
}
