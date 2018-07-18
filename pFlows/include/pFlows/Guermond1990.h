#pragma once
// Doesn't currently work:
//	At high aspect ratio, the right result is produced, but as the aspect
//  ratio decreases the reduction in lift is not as big as it should be.

#include <functional>
#include "WingProjectionGeometry.h"

namespace mFlow {
	class Guermond1990
	{
	public:
		Guermond1990();
		~Guermond1990();

		WingProjectionGeometry wing;

		// Incidence with respect to Y_Global
		std::function<double(double)> incidence;
		std::function<double(double)> zero_lift_incidence;
		// The moment about the midchord - as Eq30. This can be set
		// equal to function moment_about_midchord_no_camber_slope(double y_global)
		std::function<double(double)> moment_about_midchord;

		// The circulation at a local y coordinate on wing - Eq32
		double Gamma(double local_y);
		// The 2D lift of a the section at the local y coordinate - Eq27
		double Gamma_0(double local_y);

		// The local downwash (so not downwash in the Prandtl sense - see downwash) Eq33
		double local_downwash(double y);

		// Downwash in the Prandtl sense. Eq39 - Eq41
		double downwash(double y);

		// Compute the lift coefficient of the wing
		double lift_coefficient();

		// Possible value for moment_about_midchord:
		double moment_about_midchord_no_camber_slope(double y_global);

	private:
		// Coordinate notes: y_local is [-1, 1]

		// Eq 31 term 1, 2, 3.
		double downwash_integral1(double y);
		double downwash_integral2(double y);
		double downwash_integral3(double y);

		// Auxilery function for computing downwash integrals defined in Eq37, Eq38
		double aux_downwash_function(double y, double psi);
		double aux_downwash_function(double y);
		double daux_downwash_function_dy(double y);
	};
}
