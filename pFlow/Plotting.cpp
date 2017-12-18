#include "stdafx.h"
#include "Plotting.h"

#include "../../HBTK/HBTK/Constants.h"

void mFlow::Plotting::plot_complex_over_period(std::complex<double> plot_variable, double oscillation_period, HBTK::GnuPlot & plot, std::string line_spec)
{
	plot.hold_on();
	double omega = 2 * HBTK::Constants::pi() / oscillation_period;
	double t_max = oscillation_period;
	plot.plot([=](double t)->double {
		return (plot_variable * exp(HBTK::Constants::i() * omega * t)).real();
	}, 0., t_max, line_spec);
	return;
}
