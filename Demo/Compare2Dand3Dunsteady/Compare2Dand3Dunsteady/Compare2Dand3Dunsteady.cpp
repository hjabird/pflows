// Compare2Dand3Dunsteady.cpp : Defines the entry point for the console application.
//


#include <iostream>

#include <HBTK/Constants.h>

#include "WingProjectionGeometry.h"
#include "WingGenerators.h"
#include "Sclavounos1987.h"
#include "McGowan2011.h"
#include "Plotting.h"

int main()
{	
	std::cout << "Compare 2D / 3D unsteady sinusiodal. \n(c) HJA Bird 2018\n\n";

	double span = 1;
	double aspect_ratio = 4;
	double heave_amplitude = 0.174;
	double pitch_amplitude = 0.0;		// radians
	double phase_offset = 0.32097;		// radians
	double frequency = 103.3;		// Angular frequency (radians / s)
	double period = 2 * HBTK::Constants::pi() / frequency;
	double pitch_location = -1.0;	// LE->-1, TE->1

	
	std::complex<double> sclavounos_cl, qstat_cl, theodorsen_cl;

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, span, aspect_ratio);
	/*
	mFlow::Sclavounos1987 sclavounos;
	sclavounos.wing = wing;
	sclavounos.j = 3;
	sclavounos.U = 1;
	sclavounos.omega = frequency;
	sclavounos.number_of_terms = 8;
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
	mcgowan.pitch_location = pitch_location;
	mcgowan.semichord = wing.semichord(0);

	qstat_cl = mcgowan.qstat_unsteady();
	theodorsen_cl = mcgowan.theodorsen_unsteady();

	double span_reduced_frequency = frequency * span / 2;
	double chord_reduced_frequency = frequency * mcgowan.semichord;

	std::cout << "\nU = " << mcgowan.free_stream_vel;
	std::cout << "\nChord(y=0) = " << mcgowan.semichord * 2;
	std::cout << "\nPlunge amp = " << heave_amplitude;
	std::cout << "\nPitch amp (RAD) = " << pitch_amplitude;
	std::cout << "\nPitch amp (DEG) = " << pitch_amplitude * 180 / HBTK::Constants::pi();
	std::cout << "\nPitch phase offset (RAD) = " << phase_offset;
	std::cout << "\nPitch phase offset (DEG) = " << phase_offset * 180 / HBTK::Constants::pi();
	std::cout << "\nPitch location [-1=LE, 1=TE] = " << mcgowan.pitch_location;
	std::cout << "\nAspect ratio = " << aspect_ratio;
	std::cout << "\nPeriod = " << period;
	std::cout << "\nChord reduced freq = " << chord_reduced_frequency;
	std::cout << "\nspan reduced freq = " << span_reduced_frequency;
	std::cout << "\n\nSclavounos cl = \t" << sclavounos_cl << " -> \t" << abs(sclavounos_cl) << " at " << std::arg(sclavounos_cl);
	std::cout << "\nqstat_cl = \t\t" << qstat_cl << " -> \t" << abs(qstat_cl) << " at " << std::arg(qstat_cl);
	std::cout << "\ntheodorson cl = \t" << theodorsen_cl << " -> \t" << abs(theodorsen_cl) << " at " << std::arg(theodorsen_cl);
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

    return 0;
}

