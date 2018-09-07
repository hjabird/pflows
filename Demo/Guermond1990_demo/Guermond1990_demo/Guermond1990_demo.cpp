
#include <iostream>
#include <vector>

#include <HBTK/Constants.h>
#include <HBTK/Generators.h>
#include <HBTK/GnuPlot.h>
#include <HBTK/Paths.h>

#include "Guermond1990.h"
#include "WingGenerators.h"

int main(int argc, char* argv[]) {
	std::cout << "Demo of Guermond 1990.\n";
	std::cout << "Exe path: " << HBTK::Paths::executable_path() << "\n";
	std::cout << "Working path: " << HBTK::Paths::current_working_directory() << "\n";

	mFlow::WingProjectionGeometry wing;
	mFlow::WingGenerators::elliptic(wing, 1, 1000);
	//wing.m_TE_expr = [=](double x)->double { return wing.m_TE_expr(x) + 0.2 * x * x; };
	//wing.m_LE_expr = [=](double x)->double { return wing.m_LE_expr(x) + 0.2 * x * x; };

	HBTK::GnuPlot plt;
	wing.add_to_plot(plt);

	mFlow::Guermond1990 sim;
	sim.wing = wing;
	sim.incidence = [](double x) { return HBTK::Constants::degrees_to_radians(1.); };
	sim.zero_lift_incidence = [](double x) { return 0.0; };
	sim.moment_about_midchord = [&](double x) { return sim.moment_about_midchord_no_camber_slope(x); };

	std::cout << "Circulation at middle of wing is " << sim.Gamma_0(0.0) << "\n";
	std::cout << "Lift coeff is " << sim.lift_coefficient() << "\n";

	std::vector<double> y_values = HBTK::linspace(0, 0.99, 50);
	std::vector<double> circ(y_values.size());
	for (int i = 0; i < (int)y_values.size(); i++) {
		circ[i] = sim.Gamma(y_values[i]);
	}
	plt.hold_on();
	plt.plot(y_values, circ);

	return 0;
}