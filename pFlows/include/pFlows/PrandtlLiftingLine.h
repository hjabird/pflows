
#include <functional>

#include "WingProjectionGeometry.h"

namespace mFlow {
	class PrandtlLiftingLine {
	public:
		WingProjectionGeometry wing;

		// Incidence with respect to Y_Global
		std::function<double(double)> incidence;
		std::function<double(double)> zero_lift_incidence;

		int n_terms;

		double lift_coefficient();

	private:
		double y_to_theta(double y) const;
		double theta_to_y(double theta);

		std::vector<double> generate_collocation_points();

		double integral(int k, double theta) const;
	};
}