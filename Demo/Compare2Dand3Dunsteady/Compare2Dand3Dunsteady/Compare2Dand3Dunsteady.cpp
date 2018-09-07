// Compare2Dand3Dunsteady.cpp : Defines the entry point for the console application.
//


#include <iostream>
#include <string>

#include <HBTK/Constants.h>
#include <HBTK/CsvWriter.h>
#include <HBTK/DoubleTable.h>
#include <HBTK/Generators.h>
#include <HBTK/Integrators.h>
#include <HBTK/Paths.h>

#include "WingProjectionGeometry.h"
#include "WingGenerators.h"
#include "Sclavounos1987.h"
#include "McGowan2011.h"
#include "Plotting.h"

int main()
{
	std::cout << "Compare 2D / 3D unsteady sinusiodal. \n(c) HJA Bird 2018\n\n";
	std::cout << "Working dir: " << HBTK::Paths::current_working_directory() << "\n";
	std::cout << "Exe dir: " << HBTK::Paths::executable_path() << "\n";

	double span = 0.0762 * 12; 
	double aspect_ratio = 12;
	double heave_amplitude = 0.05 * 0.0762;// 0.05 * 0.3048 / 4;
	double pitch_amplitude = 0;		// radians
	double phase_offset = 1.5708;		// radians
	double frequency = 103.3;		// Angular frequency (radians / s)
	double period = 2 * HBTK::Constants::pi() / frequency;
	double pitch_location = 0.0;	// LE->-1, TE->1 at center span.

	std::vector<double> fq_range = HBTK::linspace(1, 20, 19);
	std::vector<double> ar_range = HBTK::linspace(1, 10, 1);
	std::string output_path = "rect_comparison.csv";

	std::complex<double> sclavounos_cl, qstat_cl, theodorsen_cl;

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, span, aspect_ratio);
	/*
	mFlow::Sclavounos1987 sclavounos;
	sclavounos.wing = wing;
	sclavounos.j = 3;
	sclavounos.U = 1;
	sclavounos.omega = frequency;
	sclavounos.number_of_terms = 16;
	sclavounos.compute_solution();
	auto pitch_complex = heave_amplitude * sclavounos.compute_equivalent_pitch_rectangular_wing(wing.semichord(0)*pitch_location);
	pitch_amplitude = - abs(pitch_complex);
	phase_offset = std::arg(pitch_complex);
	*/
	sclavounos_cl = mFlow::Sclavounos1987::conventional_lift_coefficient(heave_amplitude, pitch_amplitude, phase_offset,
		frequency, wing, mFlow::rectangular_added_mass_coefficient(span, wing.chord(0)), wing.semichord(0)*pitch_location);

	mFlow::McGowan2011 mcgowan;
	mcgowan.free_stream_vel = 1;
	mcgowan.phase_offset = phase_offset;
	mcgowan.pitch_amplitude = pitch_amplitude;
	mcgowan.plunge_amplitude = heave_amplitude;
	mcgowan.frequency = frequency;
	auto setup_2D = [&](double y){		
		double chord = wing.chord(y);
		double chord0 = wing.chord(0);
		double pitch_loc = pitch_location * chord0 / chord;
		mcgowan.pitch_location = pitch_loc;
		mcgowan.semichord = chord / 2; 
	};
	auto theodorsen_integrand = [&](double y)->std::complex<double> {
		setup_2D(y);
		return mcgowan.theodorsen_unsteady() * wing.chord(y);
	};
	auto qstat_integrand = [&](double y)->std::complex<double> {
		setup_2D(y);
		return mcgowan.qstat_unsteady() * wing.chord(y);
	};
	qstat_cl = HBTK::adaptive_simpsons_integrate(qstat_integrand, 1e-3, -wing.semispan() + 1e-5, wing.semispan() - 1e-5) / wing.area();
	theodorsen_cl = HBTK::adaptive_simpsons_integrate(theodorsen_integrand, 1e-3, -wing.semispan() + 1e-5, wing.semispan() - 1e-5) / wing.area();


	double span_reduced_frequency = frequency * span / 2;
	double chord_reduced_frequency = frequency * wing.semichord(0);

	std::cout << "\nU = " << mcgowan.free_stream_vel;
	std::cout << "\nChord(y=0) = " << wing.chord(0);
	std::cout << "\nSpan = " << wing.span;
	std::cout << "\nPlunge amp = " << heave_amplitude;
	std::cout << "\nPitch amp (RAD) = " << pitch_amplitude;
	std::cout << "\nPitch amp (DEG) = " << pitch_amplitude * 180 / HBTK::Constants::pi();
	std::cout << "\nPitch phase offset (RAD) = " << phase_offset;
	std::cout << "\nPitch phase offset (DEG) = " << phase_offset * 180 / HBTK::Constants::pi();
	std::cout << "\nPitch location [-1=LE, 1=TE] = " << mcgowan.pitch_location;
	std::cout << "\nAspect ratio = " << aspect_ratio;
	std::cout << "\nOmega = " << frequency;
	std::cout << "\nPeriod = " << period;
	std::cout << "\nChord reduced freq (y=0) = " << chord_reduced_frequency;
	std::cout << "\nSpan reduced freq = " << span_reduced_frequency;
	std::cout << "\n\nSclavounos cl = \t" << sclavounos_cl << " -> \t" << std::abs(sclavounos_cl) << " at " << std::arg(sclavounos_cl);
	std::cout << "\nqstat_cl = \t\t" << qstat_cl << " -> \t" << std::abs(qstat_cl) << " at " << std::arg(qstat_cl);
	std::cout << "\ntheodorson cl = \t" << theodorsen_cl << " -> \t" << std::abs(theodorsen_cl) << " at " << std::arg(theodorsen_cl);
	std::cout << "\n\n";

	HBTK::GnuPlot plot;
	plot.replot_off();
	plot.title("Lift against time: frequency = " + std::to_string(frequency) + ",\\n span reduced frequency = " +
		std::to_string(span_reduced_frequency) + ", chord reduced frequency = " + std::to_string(chord_reduced_frequency));
	plot.xlabel("Time (s)");
	plot.ylabel("C_l");
	plot.legend({ "QSTAT", "Theodorsen", "Sclavounos" });
	mFlow::Plotting::plot_complex_over_period(qstat_cl, period, plot, "r-");
	mFlow::Plotting::plot_complex_over_period(theodorsen_cl, period, plot, "b-");
	mFlow::Plotting::plot_complex_over_period(sclavounos_cl, period, plot, "k-");
	plot.axis_equal_off();
	plot.replot();

	if ((int)fq_range.size() > 0 && (int)ar_range.size() > 0) {
		std::cout << "Working on " << output_path << "...\n";
		std::cout << "\t" << (int)fq_range.size() << " Fq values\n";
		std::cout << "\t" << (int)ar_range.size() << " AR values\n";
		HBTK::DoubleTable table;
		table.add_column("Frequency (rad)");
		table.add_column("Aspect ratio");
		table.add_column("Chord reduced freq");
		table.add_column("Span reduced freq");
		table.add_column("Heave amp / c");
		table.add_column("Pitch amp (rad)");
		table.add_column("Phase off (rad)");
		table.add_column("Pitch loc");
		table.add_column("Theo Abs(Cl)");
		table.add_column("Theo Ph(Cl)");
		table.add_column("Scl Abs(Cl)");
		table.add_column("Scl Ph(Cl)");

		for (double fq : fq_range) {
			for (double ar : ar_range) {
				table.column(0).push_back(fq);
				mFlow::WingProjectionGeometry wing;
				mFlow::WingGenerators::rectangular(wing, span, ar);
				table.column(1).push_back(ar);
				table.column(2).push_back(fq * wing.semichord(0));
				table.column(3).push_back(fq * span / 2);
				table.column(4).push_back(heave_amplitude / wing.chord(0));
				table.column(5).push_back(pitch_amplitude);
				table.column(6).push_back(phase_offset);
				table.column(7).push_back(pitch_location);
				sclavounos_cl = mFlow::Sclavounos1987::conventional_lift_coefficient(heave_amplitude, pitch_amplitude, phase_offset,
					fq, wing, mFlow::rectangular_added_mass_coefficient(span, wing.chord(0)), wing.semichord(0)*pitch_location);
				mcgowan.frequency = fq;
				theodorsen_cl = mcgowan.theodorsen_unsteady();
				table.column(8).push_back(std::abs(theodorsen_cl));
				table.column(9).push_back(std::arg(theodorsen_cl));
				table.column(10).push_back(std::abs(sclavounos_cl));
				table.column(11).push_back(std::arg(sclavounos_cl));
			}
		}

		HBTK::CsvWriter writer;
		writer.write(std::ofstream(output_path), table);
	}
	

    return 0;
}

