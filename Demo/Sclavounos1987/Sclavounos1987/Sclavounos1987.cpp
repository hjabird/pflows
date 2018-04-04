// Sclavounos1987.cpp : Defines the entry point for the console application.
//


#include <functional>
#include <iostream>

#include <HBTK/Constants.h>
#include <HBTK/Generators.h>
#include <HBTK/GnuPlot.h>

#include "Common.h"
#include "Sclavounos1987.h"
#include "WingGenerators.h"

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
	const double max_fq_term = 16;
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
	fig7.autoscale_on();
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
	fig8.autoscale_on();
	fig8.xlabel("wd/U");
	fig8.ylabel("ph(C_L)");
	fig8.yrange(-270, 0);
	fig8.plot(fq_term, phase_res, "k-");

	return;
}

void run_elliptic_pitch_experiment() {
	mFlow::Sclavounos1987 analysis;
	mFlow::WingGenerators::elliptic(analysis.wing, 4, 4);
	analysis.U = 1;
	analysis.omega = 0.001;
	analysis.j = 5;
	analysis.number_of_terms = 8;

	std::cout << "\n\n";
	std::cout << "Computing results for comparison to Sclavounos1987 figures 7 & 8.\n";
	std::cout << "Using rectangular wing of AR = " << analysis.wing.aspect_ratio() <<
		" and span = " << analysis.wing.span << ".\n";
	std::cout << "Using " << analysis.number_of_terms << " fourier terms.\n";
	std::cout << "Using j = " << analysis.j << ".\n";
	std::cout << "Using U = " << analysis.U << ".\n\n";

	const int number_of_points = 100;
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
			this_analysis.compute_lift_coeff_j5();
		cl_real[idx] = C_L.real();
		cl_imag[idx] = C_L.imag();
		cl_abs[idx] = abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] = tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig_cl, fig_phase;
	fig_cl.title("Modulus of the pitch lift-coeffient for elliptic wing of AR = 4");
	fig_cl.xlabel("wd/U");
	fig_cl.ylabel("abs(C_L)");
	fig_cl.yrange(0, 8);
	fig_cl.hold_on();
	fig_cl.plot(fq_term, cl_abs, "k-");
	fig_cl.plot(fq_term, cl_real, "r-");
	fig_cl.plot(fq_term, cl_imag, "b-");
	fig_cl.legend({ "abs(Cl)", "Real(Cl)", "Imag(Cl)" });
	fig_cl.replot();

	fig_phase.title("Phase of pitch lift coefficeint for elliptic wing of AR = 4");
	fig_phase.xlabel("wd/U");
	fig_phase.ylabel("ph(C_L)");
	fig_phase.yrange(-270, 0);
	fig_phase.plot(fq_term, phase_res, "k-");

	system("pause");
	return;
}

void run_rectangluar_pitch_experiment() {
	mFlow::Sclavounos1987 analysis;
	mFlow::WingGenerators::rectangular(analysis.wing, 4, 4);
	analysis.U = 1;
	analysis.omega = 0.001;
	analysis.j = 5;
	analysis.number_of_terms = 8;

	std::cout << "\n\n";
	std::cout << "Computing results for comparison to Sclavounos1987 figures 7 & 8.\n";
	std::cout << "Using rectangular wing of AR = " << analysis.wing.aspect_ratio() <<
		" and span = " << analysis.wing.span << ".\n";
	std::cout << "Using " << analysis.number_of_terms << " fourier terms.\n";
	std::cout << "Using j = " << analysis.j << ".\n";
	std::cout << "Using U = " << analysis.U << ".\n\n";

	const int number_of_points = 100;
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
			this_analysis.compute_lift_coeff_j5();
		cl_real[idx] = C_L.real();
		cl_imag[idx] = C_L.imag();
		cl_abs[idx] = abs(C_L);
		phase_res[idx] = atan2(C_L.imag(), C_L.real()) / HBTK::Constants::pi() * 180;
		fq_term[idx] = tmp_omega * analysis.wing.semispan() / analysis.U;
	}

	HBTK::GnuPlot fig_cl, fig_phase;
	fig_cl.title("Modulus of the pitch lift-coeffient for rectangular wing of AR = 4");
	fig_cl.autoscale_on();
	fig_cl.xlabel("wd/U");
	fig_cl.ylabel("abs(C_L)");
	fig_cl.yrange(0, 8);
	fig_cl.hold_on();
	fig_cl.plot(fq_term, cl_abs, "k-");
	fig_cl.plot(fq_term, cl_real, "r-");
	fig_cl.plot(fq_term, cl_imag, "b-");
	fig_cl.legend({ "abs(Cl)", "Real(Cl)", "Imag(Cl)" });
	fig_cl.replot();

	fig_phase.title("Phase of pitch lift coefficeint for rectangular wing of AR = 4");
	fig_phase.autoscale_on();
	fig_phase.xlabel("wd/U");
	fig_phase.ylabel("ph(C_L)");
	fig_phase.yrange(-270, 0);
	fig_phase.plot(fq_term, phase_res, "k-");

	return;
}

void run_circulatory_lift_distribution_demo_elliptic() {
	mFlow::Sclavounos1987 analysis;
	mFlow::WingGenerators::elliptic(analysis.wing, 1, 8);
	analysis.U = 1;
	analysis.omega = 2;
	analysis.j = 3;
	analysis.number_of_terms = 8;

	const int number_of_points = 100;
	const double max_fq_term = 8;
	std::vector<double> y_pos, cl, cl_HFC, cl_HFS;
	y_pos = HBTK::linspace(-0.499, 0.499, number_of_points);
	cl.resize(number_of_points);
	cl_HFC.resize(number_of_points);
	cl_HFS.resize(number_of_points);

	auto get_Cl2D = [](double y, mFlow::Sclavounos1987 analysis) {
		std::complex<double> t1 = HBTK::Constants::pi() * (1. - analysis.F(y));
		double semichord = analysis.wing.semichord(y);
		std::complex<double> C = mFlow::Common::theodorsen_function(semichord * analysis.omega / analysis.U);
		return -4.0 * t1 * C * semichord / analysis.wing.area();
	};

	auto analysis_HFS = analysis;
	mFlow::WingGenerators::elliptic(analysis_HFS.wing, 1, 80);
	analysis_HFS.omega = 20;
	auto analysis_HFC = analysis;
	mFlow::WingGenerators::elliptic(analysis_HFC.wing, 1, 1);

	analysis.compute_solution();
	analysis_HFC.compute_solution();
	analysis_HFS.compute_solution();

	for (int i = 0; i < number_of_points; i++) {
		double y = y_pos[i];
		cl[i] = abs(get_Cl2D(y, analysis));
		cl_HFC[i] = abs(get_Cl2D(y, analysis_HFC));
		cl_HFS[i] = abs(get_Cl2D(y, analysis_HFS));
	}

	HBTK::GnuPlot fig_lift;
	fig_lift.hold_on();
	fig_lift.replot_off();
	fig_lift.title("Comparison of circulatory lift coefficient with respect to span (Elliptic wing)");
	fig_lift.xlabel("Span position");
	fig_lift.ylabel("Normalised lift per unit span");
	fig_lift.plot(y_pos, cl);
	fig_lift.plot(y_pos, cl_HFC);
	fig_lift.plot(y_pos, cl_HFS);
	fig_lift.legend({ "SRF=1, CRF=1/8", "SRF=1, CRF=1", "SRF=10, CRF=1/8" });
	fig_lift.autoscale_on();
	fig_lift.replot();

	return;
}


void run_circulatory_lift_distribution_demo_rectangular() {
	mFlow::Sclavounos1987 analysis;
	mFlow::WingGenerators::rectangular(analysis.wing, 1, 8);
	analysis.U = 1;
	analysis.omega = 2;
	analysis.j = 3;
	analysis.number_of_terms = 16;

	const int number_of_points = 100;
	std::vector<double> y_pos, cl, cl_HFC, cl_HFS;
	y_pos = HBTK::linspace(-0.499, 0.499, number_of_points);
	cl.resize(number_of_points);
	cl_HFC.resize(number_of_points);
	cl_HFS.resize(number_of_points);

	auto get_Cl2D = [](double y, mFlow::Sclavounos1987 analysis) {
		std::complex<double> t1 = HBTK::Constants::pi() * (1. - analysis.F(y));
		double semichord = analysis.wing.semichord(y);
		std::complex<double> C = mFlow::Common::theodorsen_function(semichord * analysis.omega / analysis.U);
		// return -4.0 * t1 * C * semichord / analysis.wing.area();
		return t1;
	};

	auto analysis_HFS = analysis;
	mFlow::WingGenerators::rectangular(analysis_HFS.wing, 1, 80);
	analysis_HFS.omega = 20;
	auto analysis_HFC = analysis;
	mFlow::WingGenerators::rectangular(analysis_HFC.wing, 1, 1);

	analysis.compute_solution();
	analysis_HFC.compute_solution();
	analysis_HFS.compute_solution();

	for (int i = 0; i < number_of_points; i++) {
		double y = y_pos[i];
		cl[i] = abs(get_Cl2D(y, analysis));
		cl_HFC[i] = abs(get_Cl2D(y, analysis_HFC));
		cl_HFS[i] = abs(get_Cl2D(y, analysis_HFS));
	}

	HBTK::GnuPlot fig_lift;
	fig_lift.hold_on();
	fig_lift.replot_off();
	fig_lift.title("Comparison of circulatory lift coefficient with respect to span (Rectangular wing)");
	fig_lift.xlabel("Span position");
	fig_lift.ylabel("Normalised Lift per unit span");
	fig_lift.plot(y_pos, cl);
	fig_lift.plot(y_pos, cl_HFC);
	fig_lift.plot(y_pos, cl_HFS);
	fig_lift.legend({ "SRF=1, CRF=1/8", "SRF=1, CRF=1", "SRF=10, CRF=1/8" });
	fig_lift.autoscale_on();
	fig_lift.replot();

	system("pause");
	return;
}

void run_aspect_ratio_demo_rectangular() {
	mFlow::Sclavounos1987 analysis_ar_inf;
	mFlow::WingGenerators::rectangular(analysis_ar_inf.wing, 1000, 1000);
	analysis_ar_inf.U = 1;
	analysis_ar_inf.omega = 2;
	analysis_ar_inf.j = 3;
	analysis_ar_inf.number_of_terms = 16;

	const int number_of_points = 100;
	std::vector<double> y_pos, ar_inf, ar_4, ar_2, ar_1;
	y_pos = HBTK::linspace(-0.499, 0.499, number_of_points);
	ar_inf.resize(number_of_points);
	ar_4.resize(number_of_points);
	ar_2.resize(number_of_points);
	ar_1.resize(number_of_points);

	auto get_Cl2D = [](double y, mFlow::Sclavounos1987 analysis) {
		std::complex<double> t1 = HBTK::Constants::pi() * (1. - analysis.F(y));
		double semichord = analysis.wing.semichord(y);
		std::complex<double> C = mFlow::Common::theodorsen_function(semichord * analysis.omega / analysis.U);
		// return -4.0 * t1 * C * semichord / analysis.wing.area();
		return t1;
	};

	auto analysis_ar_4 = analysis_ar_inf;
	mFlow::WingGenerators::rectangular(analysis_ar_4.wing, 4, 4);
	auto analysis_ar_2 = analysis_ar_inf;
	mFlow::WingGenerators::rectangular(analysis_ar_2.wing, 2, 2);
	auto analysis_ar_1 = analysis_ar_inf;
	mFlow::WingGenerators::rectangular(analysis_ar_1.wing, 1, 1);

	analysis_ar_inf.compute_solution();
	analysis_ar_4.compute_solution();
	analysis_ar_2.compute_solution();
	analysis_ar_1.compute_solution();

	for (int i = 0; i < number_of_points; i++) {
		double y = y_pos[i];
		ar_inf[i] = abs(get_Cl2D(y * 1000, analysis_ar_inf));
		ar_4[i] = abs(get_Cl2D(y * 4, analysis_ar_4));
		ar_2[i] = abs(get_Cl2D(y * 2, analysis_ar_2));
		ar_1[i] = abs(get_Cl2D(y * 1, analysis_ar_1));
	}

	HBTK::GnuPlot fig_lift;
	fig_lift.hold_on();
	fig_lift.replot_off();
	fig_lift.title("Lift distribution over span for different AR rectangular wings");
	fig_lift.xlabel("Span position");
	fig_lift.ylabel("Normalised Lift per unit span");
	fig_lift.plot(y_pos, ar_inf);
	fig_lift.plot(y_pos, ar_4);
	fig_lift.plot(y_pos, ar_2);
	fig_lift.plot(y_pos, ar_1);
	fig_lift.legend({ "AR_{inf}", "AR_4", "AR_2", "AR_1" });
	fig_lift.autoscale_on();
	fig_lift.replot();

	return;
}

void run_aspect_ratio_demo_elliptic() {
	mFlow::Sclavounos1987 analysis_ar_inf;
	mFlow::WingGenerators::elliptic(analysis_ar_inf.wing, 1000, 1000);
	analysis_ar_inf.U = 1;
	analysis_ar_inf.omega = 2;
	analysis_ar_inf.j = 3;
	analysis_ar_inf.number_of_terms = 16;

	const int number_of_points = 100;
	std::vector<double> y_pos, ar_inf, ar_4, ar_2, ar_1;
	y_pos = HBTK::linspace(-0.499, 0.499, number_of_points);
	ar_inf.resize(number_of_points);
	ar_4.resize(number_of_points);
	ar_2.resize(number_of_points);
	ar_1.resize(number_of_points);

	auto get_Cl2D = [](double y, mFlow::Sclavounos1987 analysis) {
		std::complex<double> t1 = HBTK::Constants::pi() * (1. - analysis.F(y));
		double semichord = analysis.wing.semichord(y);
		std::complex<double> C = mFlow::Common::theodorsen_function(semichord * analysis.omega / analysis.U);
		// return -4.0 * t1 * C * semichord / analysis.wing.area();
		return t1;
	};

	auto analysis_ar_4 = analysis_ar_inf;
	mFlow::WingGenerators::elliptic(analysis_ar_4.wing, 4, 4);
	auto analysis_ar_2 = analysis_ar_inf;
	mFlow::WingGenerators::elliptic(analysis_ar_2.wing, 2, 2);
	auto analysis_ar_1 = analysis_ar_inf;
	mFlow::WingGenerators::elliptic(analysis_ar_1.wing, 1, 1);

	analysis_ar_inf.compute_solution();
	analysis_ar_4.compute_solution();
	analysis_ar_2.compute_solution();
	analysis_ar_1.compute_solution();

	for (int i = 0; i < number_of_points; i++) {
		double y = y_pos[i];
		ar_inf[i] = abs(get_Cl2D(y * 1000, analysis_ar_inf));
		ar_4[i] = abs(get_Cl2D(y * 4, analysis_ar_4));
		ar_2[i] = abs(get_Cl2D(y * 2, analysis_ar_2));
		ar_1[i] = abs(get_Cl2D(y * 1, analysis_ar_1));
	}

	HBTK::GnuPlot fig_lift;
	fig_lift.hold_on();
	fig_lift.replot_off();
	fig_lift.title("Lift distribution over span for different AR elliptic wings");
	fig_lift.xlabel("Span position");
	fig_lift.ylabel("Normalised Lift per unit span");
	fig_lift.plot(y_pos, ar_inf);
	fig_lift.plot(y_pos, ar_4);
	fig_lift.plot(y_pos, ar_2);
	fig_lift.plot(y_pos, ar_1);
	fig_lift.legend({ "AR_{inf}", "AR_4", "AR_2", "AR_1" });
	fig_lift.autoscale_on();
	fig_lift.replot();

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

	//run_elliptic_experiment();
	run_rectangular_experiment();
	//run_elliptic_pitch_experiment();
	//run_rectangluar_pitch_experiment();
	//run_circulatory_lift_distribution_demo_elliptic();
	//run_circulatory_lift_distribution_demo_rectangular();
	//run_aspect_ratio_demo_elliptic();

	return 0;
}
