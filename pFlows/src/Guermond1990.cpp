#include "Guermond1990.h"

#include <cassert>
#include <HBTK/NumericalDifferentiation.h>
#include <HBTK/GaussLegendre.h>
#include <HBTK/Integrators.h>
#include <HBTK/Constants.h>
#include <HBTK/Checks.h>


mFlow::Guermond1990::Guermond1990()
{
}


mFlow::Guermond1990::~Guermond1990()
{
}


double mFlow::Guermond1990::Gamma(double local_y)
{
	assert(abs(local_y) <= 1.);
	double global_y = local_y  * wing.semispan();
	return HBTK::Constants::pi() * wing.chord(global_y) * cos(wing.midchord_sweep_angle(global_y)) *
		(incidence(global_y) - zero_lift_incidence(global_y) + local_downwash(local_y));
}

double mFlow::Guermond1990::Gamma_0(double local_y)
{
	assert(abs(local_y) <= 1);
	double y_global = local_y * wing.semispan();
	double cos_lambda = cos(wing.midchord_sweep_angle(y_global));
	double equivalent_incidence = incidence(y_global) - zero_lift_incidence(y_global);
	assert(HBTK::check_finite(equivalent_incidence));
	return HBTK::Constants::pi() * wing.chord(y_global) * cos_lambda * equivalent_incidence;
}

double mFlow::Guermond1990::local_downwash(double y)
{
	assert(abs(y) <= 1.0);
	double A = wing.aspect_ratio();
	double c = wing.chord(y * wing.semispan());
	double K = 0.5;
	double Lambda = wing.midchord_sweep_angle(y * wing.semispan());
	double r = wing.midchord_radius_of_curvature(y * wing.semispan());
	double pi = HBTK::Constants::pi();
	double term_11 = (Gamma_0(y) / (4 * pi*r) + sin(Lambda)*
		HBTK::central_difference_O1A2([&](double y) { return Gamma_0(y); }, y) / (2 * pi)) *
		(log(c / A) - 2 * log(2) + 1 - 2 * K);
	double term_12 = (-2 / c) * (moment_about_midchord(y) / (4 * pi * r) + sin(Lambda) *
		HBTK::central_difference_O1A2([&](double y) { return moment_about_midchord(y); }, y) / (2 * pi));
	double term_13 = downwash(y);
	double term_14 = (Gamma_0(y) / (4 * pi * r)) * (1 - pow(tan(Lambda), 2) - log(2 / pow(cos(Lambda), 2)));
	double term_15 = (HBTK::central_difference_O1A2([&](double y) { return Gamma_0(y); }, y) / (2 * pi))
		* (log((1 + sin(Lambda)) / cos(Lambda)) - sin(Lambda) * log(2 / (pow(cos(Lambda), 2))));

	double return_value = (1./A) * (term_11 + term_12 + term_13 + term_14 + term_15);
	return return_value;
}


double mFlow::Guermond1990::downwash(double y)
{
	assert(abs(y) <= 1.);
	double term_1 = downwash_integral1(y) / (4. * HBTK::Constants::pi());
	double term_2 = Gamma_0(y) * downwash_integral2(y) / (4. * HBTK::Constants::pi());
	double term_3 = HBTK::central_difference_O1A2([&](double y) {return Gamma_0(y); }, y) * downwash_integral3(y)
		/ (4. * HBTK::Constants::pi());
	return term_1 + term_2 + term_3;
}

double mFlow::Guermond1990::lift_coefficient()
{
	auto integrand = [&](double y)->double {
		return Gamma(y) * 2.0;
	};
	auto quad = HBTK::gauss_legendre(50);
	double integral = quad.integrate(integrand) * wing.semispan();
	integral /= wing.area();
	return integral;
}

double mFlow::Guermond1990::moment_about_midchord_no_camber_slope(double y_global)
{
	double y = y_global / wing.span;
	return Gamma_0(y) * wing.chord(y_global) * ( 0.25 - 0.5 );
}


double mFlow::Guermond1990::downwash_integral1(double y)
{
	// First term of Eq39.
	assert(abs(y) <= 1.);
	double integral;

	double gamma0y = Gamma_0(y);
	double dgamma0y_dy = HBTK::central_difference_O1A2([&](double x) { return Gamma_0(x); }, y);

	auto integrand = [&](double psi) {
		double term_1 = (Gamma_0(psi) - gamma0y - (psi - y) * dgamma0y_dy) / pow(y - psi, 2);
		double term_2 = 1. - (psi > y ? 1. : -1.) * aux_downwash_function(y, psi);
		return term_1 * term_2;
	};
	auto quad = HBTK::gauss_legendre(20);
	integral = wing.semispan() * quad.integrate(integrand);
	assert(HBTK::check_finite(integral));
	return integral;
}

double mFlow::Guermond1990::downwash_integral2(double y)
{
	assert(abs(y) <= 1.);
	// Eq 40 with smooth wing assumptions.
	double term_1, term_2, term_3, term_4,
		term_21, term_22;

	term_1 = -2. / (1. - y * y);

	term_21 = 2. * y / (1. - y * y);
	term_22 = sin(wing.midchord_sweep_angle(y * wing.semispan()));
	term_2 = term_21 * term_22;

	term_3 = - log(1. - y * y) / (2. * wing.midchord_radius_of_curvature(wing.semispan() * y));

	auto quad = HBTK::gauss_legendre(20);
	double static_aux_term = aux_downwash_function(y);
	double static_aux_derivative = daux_downwash_function_dy(y);
	auto integrand = [&](double psi) {
		return aux_downwash_function(y, psi) - static_aux_term - (psi - y) * static_aux_derivative
			/ ((psi - y) * abs(psi - y));
	};
	term_4 = -1. * quad.integrate(integrand);
	assert(HBTK::check_finite(term_4));
	return term_1 + term_2 + term_3 + term_4;
}

double mFlow::Guermond1990::downwash_integral3(double y)
{
	// Eq 41 with smooth wing assumptions.
	assert(abs(y) <= 1.);
	double term_1, term_2, term_3;

	term_1 = log((1. - y) / (1. + y));
	term_2 = -1 * sin(wing.midchord_sweep_angle(y * wing.semispan()))
		* log(1 - y * y);

	auto quadrature = HBTK::gauss_legendre(20);
	double precomputed_downwash_aux_at_y = aux_downwash_function(y);
	auto integrand = [&](double psi) {
		return (aux_downwash_function(y, psi) - precomputed_downwash_aux_at_y) / abs(psi - y);
	};
	term_3 = -1 * quadrature.integrate(integrand);

	return term_1 + term_2 + term_3;
}


double mFlow::Guermond1990::aux_downwash_function(double y, double psi)
{
	assert(abs(y) <= 1.0);
	assert(abs(psi) <= 1.0);
	if (psi == y) return aux_downwash_function(y);
	double x0y = wing.midchord_X(y * wing.semispan()) / wing.semispan();
	double x0psi = wing.midchord_X(psi * wing.semispan()) / wing.semispan();

	double num1 = psi > y ? 1. : -1.;
	double num2 = x0psi - x0y;
	double num = num1 * num2;

	double den1 = pow(psi - y, 2);
	double den2 = pow(x0psi - x0y, 2);
	double den = sqrt(den1 + den2);

	return num / den;
}

double mFlow::Guermond1990::aux_downwash_function(double y)
{
	assert(abs(y) <= 1.);
	return sin( wing.midchord_sweep_angle(y * wing.semispan()) );
}


double mFlow::Guermond1990::daux_downwash_function_dy(double y)
{
	assert(abs(y) <= 1.);
	return 1./ (2. * wing.midchord_radius_of_curvature(y* wing.semispan()));
}
