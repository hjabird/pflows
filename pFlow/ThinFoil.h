#pragma once
/*////////////////////////////////////////////////////////////////////////////
ThinFoil.h

An object and (mostly) methods to represent a thin aerofoil in steady
flow.

Copyright 2017 HJA Bird

mFlow is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mFlow is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mFlow.  If not, see <http://www.gnu.org/licenses/>.
*/////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cassert>

#include "../../HBTK/HBTK/Integrators.h"
#include "../../HBTK/HBTK/Constants.h"

template <typename TCoeffStore>
class ThinFoil 
{
public:
	// Suggested reference is Katz and Plotkin 2001, Chapter 5.

	// Physical description coefficients.
	TCoeffStore m_coefficients;

	// The function associated with a single term - ie for n > 0: sin(n theta)
	constexpr double series_term_func(double theta, int n);
	// Evaluate a term in the series including its coeff: A_n sin(n theta) for n > 0.
	const double series_term_eval(double theta, int n);
	// Sum the entire series.
	const double series_sum(double theta);
	// Sum the entire series with additional angle of attack term.
	const double series_sum(double theta, double AoA);
	// The vorticity function as function of theta, U_infty
	const double vorticity(double theta, double U_infty, double AoA);

	// Integral function for a single term of the series without coefficient.
	// DO NOT USE except for n=0, n=1 or you'll get garbage back.
	constexpr double series_term_func_intx(double x, int n);

	// The x dir velocity in the limit of z=0+  and z=0- respectively.
	const double u_zero_plus(double x, double U_infty, double AoA);
	const double u_zero_minus(double x, double U_infty, double AoA);
	// Delta u accross foil:
	const double u_zero_delta(double x, double U_infty, double AoA);


	// The differential with respect to x of a single function: ie for n > 0 (d/dx) sin(n theta)
	constexpr double series_term_func_dx(double theta, int n);
	// Evaluate the derivative with respect to x of a term inc its coeff.
	constexpr double series_term_eval_dx(double theta, int n);
	// Sum the derivatives of the entire series
	const double series_sum_dx(double theta);
	// Sum the derivates with an additional angle of attack term.
	const double series_sum_dx(double theta, double AoA);

	// Compute the coefficients for a thin foil from a preexisting foil. The function
	// must be of the form double(double x) for x in [0, chord].
	template<typename TFunc>
	void generate_coefficients_from_foil_slope(TFunc & dydx_foil_at_x, double chord);

	// Mappings of x in [0, 1] to theta in [0, pi] and reverse.
	constexpr double theta_to_chord_pos(double theta);
	constexpr double chord_pos_to_theta(double x);
};


// DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::series_term_func(double theta, int n)
{
	assert(n >= 0);
	assert(theta >= 0);
	assert(theta <= Constants::pi() * 1.01);
	if (n = 0) {
		return (1 - cos(theta)) / sin(theta);
	}
	else {
		return (sin(theta));
	}
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::series_term_eval(double theta, int n)
{
	assert(n > 0);
	assert(n < (int)m_coefficients.size());
	return m_coefficients[n] * series_term_func(theta, n);
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::series_sum(double theta)
{
	assert(theta >= 0);
	assert(theta <= Constants::pi() * 1.01);
	double sum = 0;
	for (auto n = 0; n < m_coefficents.size(); n++) {
		sum += series_term_eval(theta, n);
	}
	return sum;
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::series_sum(double theta, double AoA)
{
	assert(theta >= 0);
	assert(theta <= Constants::pi() * 1.01);
	return series_sum(theta) - AoA * series_term_function(theta, 0);
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::vorticity(double theta, double U_infty, double AoA)
{
	return 2 * U_infty * series_sum(theta, AoA);
}

template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::series_term_func_intx(double x, int n)
{
	assert(theta >= 0);
	assert(theta <= Constants::pi());
	assert(n >= 0);

	double result;
	if (n == 0) {
		result = x - sqrt(4 * x - 4 * x*x) / 2;
	}
	else if (n == 1) {
		double tInner = 2 * x - 1;
		result = (asin(tInner) + tInner*sqrt(1 - (tInner * tInner))) / 4;
	}
	else {
		// Nasty and complicated - definite integration numerically seems like
		// an easier option.
		result = NAN;
	}

	return result;
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::u_zero_plus(double x, double U_infty, double AoA)
{
	return vorticity(chord_pos_to_theta(x), U_infty, AoA) / 2; // Katz Eq(5.37)
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::u_zero_minus(double x, double U_infty, double AoA)
{
	return - vorticity(chord_pos_to_theta(x), U_infty, AoA) / 2; // Katz Eq(5.37)
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::u_zero_delta(double x, double U_infty, double AoA)
{
	return vorticity(chord_pos_to_theta(x), U_infty, AoA); // Katz Eq(5.37)
}


template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::series_term_func_dx(double theta, int n)
{
	assert(n >= 0);
	assert(theta >= 0);
	assert(theta < Constants::pi() * 1.01);
	double x = (1 - cos(theta)) / 2;
	if (n > 0) {
		return n * cos(n*theta) / sqrt(x - x*x);
	}
	else {
		return ((1 - cos(theta)) / (sin(theta) * sin(theta))) / sqrt(x - x*x);
	}
}

template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::series_term_eval_dx(double theta, int n)
{
	assert(n >=  0);
	assert(n < (int)m_coefficients.size());
	assert(theta >= 0);
	assert(theta < Constants::pi() * 1.01);
	return m_coefficients[n] * series_term_func_dx(theta, n);
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::series_sum_dx(double theta)
{
	double sum = 0;
	for (auto idx = 0; idx < (int)m_coefficients.size(); idx++) {
		sum += series_term_eval_dx(theta, idx);
	}
	return sum;
}

template<typename TCoeffStore>
inline const double ThinFoil<TCoeffStore>::series_sum_dx(double theta, double AoA)
{
	return series_sum_dx(theta) - AoA * series_term_func_dx(theta, 0);
}

template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::theta_to_chord_pos(double theta)
{
	assert(theta >= 0);
	assert(theta <= Constants::pi() * 1.01);
	return (1 - cos(theta)) / 2.0;
}

template<typename TCoeffStore>
inline constexpr double ThinFoil<TCoeffStore>::chord_pos_to_theta(double x)
{
	assert(x >= 0);
	assert(x <= 1);
	return acos(1 - 2 * x);
}

template<typename TCoeffStore>
template<typename TFunc>
inline void ThinFoil<TCoeffStore>::generate_coefficients_from_foil_slope(TFunc & dydx_foil_at_x, double chord)
{
	assert(chord > 0);
	assert((int)m_coefficients.size() > 0);
	const int n_terms = (int)m_coefficients.size();
	const double pi = Constants::pi();

	auto integrate = [&](int term_idx)->double {
		auto integrand = [&](double theta)->double {
			return dydx_foil_at_x(theta_to_chord_pos(theta)*chord) * cos(term_idx * theta);
		};
		return Quad::adaptive_simpsons_integrate(integrand, 1e-9, 0.0, pi);
	};

	m_coefficients[0] = integrate(0)/pi;
	for (auto idx = 1; idx < (int)n_terms; idx++) {
		m_coefficients[idx] = 2 * integrate(idx) / pi;
	}

	return;
}
