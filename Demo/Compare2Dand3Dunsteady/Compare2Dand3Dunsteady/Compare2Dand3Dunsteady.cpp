// Compare2Dand3Dunsteady.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <HBTK/Constants.h>

#include "../../../pFlow/WingProjectionGeometry.h"
#include "../../../pFlow/WingGenerators.h"
#include "../../../pFlow/Sclavounos1987.h"
#include "../../../pFlow/McGowan2011.h"
#include "../../../pFlow/Plotting.h"

int main()
{	
	double span = 4;
	double aspect_ratio = 4;
	double heave_amplitude = 0.05;	
	double pitch_amplitude = 0;		// radians
	double phase_offset = 0.249;	// radians
	double frequency = 7.86;		// Angular frequency (radians / s)
	double period = 2 * HBTK::Constants::pi() / frequency;

	
	std::complex<double> sclavounos_cl, qstat_cl, theodorsen_cl;

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::rectangular(wing, span, aspect_ratio);

	mFlow::Sclavounos1987 sclavounos;
	sclavounos.wing = wing;
	sclavounos.j = 3;
	sclavounos.U = 1;
	sclavounos.omega = frequency;
	sclavounos.number_of_terms = 8;
	sclavounos.compute_solution();
	auto pitch_complex = heave_amplitude * sclavounos.compute_equivalent_pitch_rectangular_wing(-0.0);
	pitch_amplitude = -abs(pitch_complex);
	phase_offset = std::arg(pitch_complex);

	sclavounos_cl = mFlow::Sclavounos1987::conventional_lift_coefficient(heave_amplitude, pitch_amplitude, phase_offset,
		frequency, wing, mFlow::rectangular_added_mass_coefficient(span, wing.chord(0)));

	mFlow::McGowan2011 mcgowan;
	mcgowan.free_stream_vel = 1;
	mcgowan.phase_offset = phase_offset;
	mcgowan.pitch_amplitude = pitch_amplitude;
	mcgowan.plunge_amplitude = heave_amplitude;
	mcgowan.frequency = frequency;
	mcgowan.pitch_location = 0;
	mcgowan.semichord = wing.semichord(0);

	qstat_cl = mcgowan.qstat_unsteady();
	theodorsen_cl = mcgowan.theodorsen_unsteady();

	double span_reduced_frequency = frequency * span / 2;
	double chord_reduced_frequency = frequency * mcgowan.semichord;

	HBTK::GnuPlot plot;
	plot.replot_off();
	plot.title("Lift against time: frequency = " + std::to_string(frequency) + ",\\n span_reduced_frequency = " +
		std::to_string(span_reduced_frequency) + ", chord_reduced_frequency = " + std::to_string(chord_reduced_frequency));
	plot.xlabel("Time (s)");
	plot.ylabel("C_l");
	plot.legend({ "QSTAT", "Theodorsen", "Sclavounos" });
	mFlow::Plotting::plot_complex_over_period(qstat_cl, period, plot, "r-");
	mFlow::Plotting::plot_complex_over_period(theodorsen_cl, period, plot, "b-");
	mFlow::Plotting::plot_complex_over_period(sclavounos_cl, period, plot, "k-");
	plot.replot();

    return 0;
}

