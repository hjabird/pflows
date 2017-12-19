#include "stdafx.h"
#include "Sclavounos1987.h"
/*////////////////////////////////////////////////////////////////////////////
Sclavounos1987.cpp

Equations from the paper "An unsteady lifting-line theory", P.D. Sclavounos,
1987, Journal of Engineering Mathmatics.
https://doi.org/10.1007/BF00127464

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

#include <array>
#include <algorithm>
#include <numeric>
#include <functional>
#include <tuple>

#include "../../HBTK/HBTK/Constants.h"
#include "../../HBTK/HBTK/GaussLegendre.h"
#include "../../HBTK/HBTK/Remaps.h"
#include "../../HBTK/HBTK/Integrators.h"
#include "../../HBTK/HBTK/Checks.h"
#include "../../HBTK/HBTK/Generators.h"
#include "../../HBTK/HBTK/GnuPlot.h"
#include "../../HBTK/HBTK/Interpolators.h"

#include <Eigen/LU>

#include "Common.h"

namespace mFlow {
	Sclavounos1987::Sclavounos1987()
	{
	}


	Sclavounos1987::~Sclavounos1987()
	{
	}

	std::complex<double> Sclavounos1987::conventional_lift_coefficient(double heave_amplitude, 
		double pitch_amplitude, double phase_offset, double frequency, WingProjectionGeometry wing,
		double heave_added_mass_a33)
	{
		assert(frequency > 0);
		assert(HBTK::check_finite(heave_amplitude));
		assert(HBTK::check_finite(pitch_amplitude));

		std::complex<double> heave_cl, pitch_cl;
		Sclavounos1987 analysis;
		analysis.wing = wing;
		analysis.U = 1.0;
		analysis.omega = frequency;
		analysis.number_of_terms = 8;

		analysis.j = 3;
		analysis.compute_solution();
		heave_cl = analysis.compute_lift_coeff_j3(heave_added_mass_a33) * heave_amplitude  
			* (analysis.omega / analysis.U);

		analysis.j = 5;
		analysis.compute_solution();
		pitch_cl = analysis.compute_lift_coeff_j5() * pitch_amplitude * (analysis.omega/analysis.U);
		pitch_cl *= exp(HBTK::Constants::i() * phase_offset);

		assert(HBTK::check_finite(pitch_cl));
		assert(HBTK::check_finite(heave_cl));
		return (heave_cl + pitch_cl) * HBTK::Constants::i();
	}


	void Sclavounos1987::set_collocation_points()
	{
		// After Eq 7.1 states these as local maxima of sin((N+1)theta). This doesn't work.
		m_collocation_points.resize(number_of_terms);
		const auto hpi = HBTK::Constants::pi() / 2;
		for (auto idx = 0; idx < number_of_terms; idx++) {
			m_collocation_points[idx] = (hpi + idx * HBTK::Constants::pi()) / (2 * number_of_terms);
			assert(m_collocation_points[idx] <= HBTK::Constants::pi() && m_collocation_points[idx] >= 0.);
		}
		return;
	}

	Eigen::Matrix<std::complex<double>, -1, 1> Sclavounos1987::rhs_vector_of_circulations()
	{
		assert((j == 3) || (j == 5));

		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> circ_vector;
		circ_vector.resize(number_of_terms);
		for (int i = 0; i < number_of_terms; i++) {
			auto y = wing.semispan()*cos(m_collocation_points[i]);
			if (j == 3) {
				circ_vector[i] = d_3(y);
			}
			else if (j == 5) { 
				circ_vector[i] = d_5(y) - U * d_3(y) / (HBTK::Constants::i() * omega); 
			}
		}
		return circ_vector;
	}

	Eigen::MatrixXd Sclavounos1987::gamma_terms_matrix()
	{
		assert((int)m_collocation_points.size() == number_of_terms);
		Eigen::MatrixXd gamma_matrix; // The matrix that directly gives vorticities at given points.
		gamma_matrix.resize(number_of_terms, number_of_terms);

		for (int i = 0; i < number_of_terms; i++) {
			for (int j = 0; j < number_of_terms; j++) {
				assert(m_collocation_points[i] < HBTK::Constants::pi());
				assert(m_collocation_points[i] > 0);
				gamma_matrix(i, j) = sin((2 * j + 1) * m_collocation_points[i]);
			}
		}
		return gamma_matrix;
	}

	std::complex<double> Sclavounos1987::integral_diff_coefficient_for_finding_circ(double y_position)
	{
		std::complex<double> ext_coeff;
		if (omega != 0) {
			ext_coeff = d_3(y_position) / (2 * HBTK::Constants::pi() * omega * HBTK::Constants::i());
		}
		else {
			ext_coeff = -wing.semichord(y_position);
		}
		return ext_coeff;
	}


	std::complex<double> Sclavounos1987::d_3(double y)
	{
		assert(y <= abs(wing.semispan()));
		auto l = wing.semichord(y);		// semichord
		auto v = omega / U;				// Normalised frequency 

		if (omega == 0.0) { return -2 * HBTK::Constants::pi() * l; } // See Eq(5.5)

		auto numerator = 4. * U * exp( - HBTK::Constants::i() * l * v);
		auto denominator = HBTK::Constants::i() * Common::hankel2_0(v * l) + Common::hankel2_1(v * l);
		assert(HBTK::check_finite(numerator));
		assert(HBTK::check_finite(denominator));
		return numerator / denominator;
	}


	std::complex<double> Sclavounos1987::d_5(double y)
	{
		// Eq4.8. Input variable y is checked in d_3 function.
		return - 0.5 * wing.semichord(y) * d_3(y); 
	}


	std::complex<double> Sclavounos1987::K(double y)
	{
		assert(y != 0.);

		if (omega == 0.) { return K_term1(y); } // Eq 5.4
		else { return K_term1(y) + K_term2(y) + K_term3(y); } // Full K.
	}


	double Sclavounos1987::K_term1(double y)
	{
		assert(y != 0);
		return K_term1_numerator(y) * K_term1_singularity(y);
	}


	double Sclavounos1987::K_term1_singularity(double y)
	{
		assert(y != 0);
		return 1. / y;
	}


	double Sclavounos1987::K_term1_numerator(double y)
	{
		return 0.5 * exp(- (omega / U) * abs(y));
	}


	std::complex<double> Sclavounos1987::K_term2(double y)
	{
		// Avoid evaluation only to multiply by zero.
		if (omega == 0) { return 0; }

		auto coeff = 0.5 * (y >= 0 ? 1. : -1.);
		auto nu = omega / U;
		auto E_1_term = -HBTK::Constants::i() * nu 
			* Common::exponential_int_E1(nu * abs(y));

		assert(HBTK::check_finite(E_1_term));
		return coeff * -1. * E_1_term;
	}


	std::complex<double> Sclavounos1987::K_term3(double y)
	{
		// Avoid evaluation only to multiply by zero.
		if (omega == 0) { return 0; } 

		auto coeff = -0.5 * (y >= 0 ? 1. : -1.);
		auto nu = omega / U;
		auto P_term = nu * P(nu * abs(y));
		
		assert(HBTK::check_finite(P_term));
		return coeff * P_term;
	}


	std::complex<double> Sclavounos1987::integrate_gammaprime_K(double y_position, int k)
	{
		assert(abs(y_position) <= wing.semispan());
		assert(k >= 0);

		// We're evaluating int( Gamma'(eta) * K(y-eta deta)
		// K can be expanded into 3 terms as considered separatively as follows:
		
		auto integral1 = integrate_gammaprime_K_term1(y_position, k);
		auto integral2 = integrate_gammaprime_K_term2(y_position, k);
		auto integral3 = integrate_gammaprime_K_term3(y_position, k);

		return integral1 + integral2 + integral3;
	}


	std::complex<double> Sclavounos1987::integrate_gammaprime_K_term1(double y_position, int k)
	{
		auto theta_sing = acos(y_position / wing.semispan());
		
		const int num_quad_points = 30;
		auto quad = split_quad(num_quad_points, theta_sing);

		// We're using the singularity subtraction method.
		auto singularity_coefficient = dsintheta_dtheta(theta_sing, k) * K_term1_numerator(0);

		auto numerical_integrand = [&](double theta0){
			auto eta = wing.semispan() * cos(theta0);
			auto singular_part = K_term1_singularity(y_position - eta);
			auto nonsingular_K = K_term1_numerator(y_position - eta);
			auto gamma_dtheta = dsintheta_dtheta(theta0, k);

			auto singularity_subtraction = 
				(nonsingular_K * gamma_dtheta - singularity_coefficient);
			return singular_part * singularity_subtraction;
		};

		auto integral = HBTK::static_integrate(numerical_integrand, 
					std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size())
			+ HBTK::static_integrate(numerical_integrand, 
					std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size())
			+ singularity_coefficient * 0.; // Glauert integral

		assert(HBTK::check_finite(integral));
		// We switched the limits to integrate 0 -> pi, so we need a negative sign here.
		return -integral;
	}


	std::complex<double> Sclavounos1987::integrate_gammaprime_K_term2(double y, int k)
	{
		assert(abs(y) <= wing.semispan());
		assert(k >= 0);

		// The singularity position in terms of the angle theta.
		auto theta_sing = acos(y / wing.semispan());
		double v = omega / U;

		auto integral_coefficient = HBTK::Constants::i() * omega / (2 * U);

		// The non-singular part of the integrand
		auto non_singular = [&](double theta)->double {
			return (2 * k + 1) * cos((2 * k + 1)*theta) / (v * wing.semispan() * sin(theta));
		};
		// The singular part of the integrand
		auto singular = [&](double theta)->std::complex<double> {
			return v * wing.semispan() * sin(theta) * (theta > theta_sing ? 1. : -1.)
				* Common::exponential_int_E1(v * wing.semispan() * abs(cos(theta_sing) - cos(theta)));
		};
		// The singular part of the integrand evaluated.
		auto singular_integral = v * wing.semispan() * (
			(cos(theta_sing)+1.) * Common::exponential_int_E1(v * wing.semispan() * (cos(theta_sing)+1.))
			+ (cos(theta_sing) -1.) * Common::exponential_int_E1(v * wing.semispan() * (1. - cos(theta_sing)))) 
			+ exp(v * wing.semispan()*(cos(theta_sing) - 1.)) - exp(-v * wing.semispan()*(cos(theta_sing) + 1.));

		// The singularity subtraction method is used here.
		auto ssm_variable = non_singular(theta_sing);
		auto numerical_integrand = [&](double theta)->std::complex<double> {
			auto singular_var = singular(theta);
			auto non_singular_var = non_singular(theta);
			auto singularity_subtraction = non_singular_var - ssm_variable;
			auto integrand = singular_var * singularity_subtraction;
			return integrand;
		};

		auto quad = split_quad(50, theta_sing); // points_lower, weights_lower, points_upper, weights_upper
		auto int_lower = HBTK::static_integrate(numerical_integrand, 
			std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size());
		auto int_upper = HBTK::static_integrate(numerical_integrand, 
			std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size());

		auto complete_integral = integral_coefficient * (int_lower + int_upper + ssm_variable * singular_integral);
		assert(HBTK::check_finite(complete_integral));
		return complete_integral;
	}


	std::complex<double> Sclavounos1987::integrate_gammaprime_K_term3(double y, int k)
	{
		// Integrating the P term of of int ( Gamma'(eta) K(y-eta) deta).
		// We're integrating in [0, singularity], [singularity, pi]
		auto theta_sing = acos(y / wing.semispan());
		auto quad = split_quad(40, theta_sing); // lower_points, lower_weights, upper_points, upper_weights
		auto integrand = [&](double theta0)->std::complex<double> {
			auto eta = wing.semispan() * cos(theta0);
			return dsintheta_dtheta(theta0, k) * K_term3(y - eta);
		};

		auto integral = HBTK::static_integrate(integrand, std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size())
			+ HBTK::static_integrate(integrand, std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size());

		assert(HBTK::check_finite(integral));
		return integral;
	}


	std::complex<double> Sclavounos1987::F(double y)
	{
		assert(y <= abs(wing.semispan()));
		assert((j == 3) || (j == 5));
		if (j == 3) {
			return 1. - solution_vorticity(y) / d_3(y);
		}
		else { //j == 5
			return (d_5(y) - solution_vorticity(y)) / d_3(y) - U / (HBTK::Constants::i() * omega);
		}
	}


	double Sclavounos1987::dtheta_dy(double y)
	{
		assert( abs(y) <= wing.semispan() );
		auto result = -1. / sqrt(pow(wing.semispan(), 2) - pow(y, 2));
		assert(HBTK::check_finite(result));
		return result;
	}


	double Sclavounos1987::dsintheta_dy(double y, int k)
	{
		assert(y <= abs(wing.semispan()));
		assert(k >= 0);
		
		auto theta = acos(y / wing.semispan());
		auto dtdy = dtheta_dy(y);
		auto dGammadt = dsintheta_dtheta(theta, k);

		return dGammadt * dtdy;
	}


	double Sclavounos1987::dsintheta_dtheta(double theta, int k)
	{
		assert(theta >= 0);
		assert(theta <= HBTK::Constants::pi());
		auto dGamma_dt = (2 * k + 1) * cos((2 * k + 1) * theta);
		return dGamma_dt;
	}


	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
		Sclavounos1987::split_quad(int total_pts, double split_theta)
	{
		assert(total_pts >= 2);

		// Compute the number of points for each quadrature based on linear split.
		int n_pts_lower = (int)ceil((total_pts-2) * split_theta / HBTK::Constants::pi()) + 1;
		int n_pts_upper = (total_pts-1) - n_pts_lower + 1;

		std::vector<double> points_upper, weights_upper, points_lower, weights_lower;
		points_upper.resize(n_pts_upper);	weights_upper.resize(n_pts_upper);
		points_lower.resize(n_pts_lower);	weights_lower.resize(n_pts_lower);
		HBTK::gauss_legendre<double>(n_pts_lower, points_lower, weights_lower);
		HBTK::gauss_legendre<double>(n_pts_upper, points_upper, weights_upper);
		for (int idx = 0; idx < (int)n_pts_lower; idx++) {
			HBTK::linear_remap(points_lower[idx], weights_lower[idx], -1., 1., 0., split_theta);
		}
		for (int idx = 0; idx < (int)n_pts_upper; idx++) {
			HBTK::linear_remap(points_upper[idx], weights_upper[idx], -1., 1., split_theta, HBTK::Constants::pi());
		}

		return std::make_tuple(points_lower, weights_lower, points_upper, weights_upper);
	}


	void Sclavounos1987::compute_solution()
	{
		assert(U > 0);
		assert(omega >= 0);
		assert(j == 3 || j == 5);
		assert(number_of_terms > 0);

		set_collocation_points();
		m_solution.resize(number_of_terms);
		auto gamma_matrix = gamma_terms_matrix();

		// Compute integ_diff_matrix - the second term of eq5.3 with the integral in it.
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> integ_diff_matrix;
		integ_diff_matrix.resize(number_of_terms, number_of_terms);
		for (int i = 0; i < number_of_terms; i++) {
			auto y_position = wing.semispan() * cos(m_collocation_points[i]);
			std::complex<double> ext_coeff = integral_diff_coefficient_for_finding_circ(y_position);
			for (int j = 0; j < number_of_terms; j++) 
			{ 
				integ_diff_matrix(i, j) = ext_coeff * integrate_gammaprime_K(y_position, j);
			} // End For in j
		} // End for in i

		auto RHS_vector = Sclavounos1987::rhs_vector_of_circulations();
		m_solution = (gamma_matrix - integ_diff_matrix).lu().solve(RHS_vector);
		return;
	}


	std::complex<double> Sclavounos1987::solution_vorticity(double y)
	{
		assert(m_solution.size() == number_of_terms);
		assert(y <= abs(wing.semispan()) );

		std::complex<double> sum(0, 0);
		double theta;
		for (int idx = 0; idx < number_of_terms; idx++) {
			theta = acos(y / wing.semispan());
			assert(HBTK::check_finite(m_solution(idx)));
			sum += sin((2 * idx + 1) * theta) * m_solution(idx);
		}

		assert(HBTK::check_finite(sum));
		return sum;
	}


	std::complex<double> Sclavounos1987::compute_lift_coeff_j3(double heave_added_mass)
	{
		assert(omega != 0.);
		assert(j == 3);
		std::complex<double> term_1, term_2,
			term_11, term_12, term_21, term_22;

		auto integrand = [&](double theta) -> std::complex<double> {
			auto y = wing.semispan() * cos(theta);
			auto semichord = wing.semichord(y);
			auto F_3 = F(y);
			auto C = Common::theodorsen_function((omega / U) * semichord);
			return C * semichord * HBTK::Constants::pi() * (1. - F_3) * sin(theta);
		};

		term_11 = -8. * wing.semispan() / wing.area();
		const int n_pts = 40;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			HBTK::linear_remap(points[idx], weights[idx], -1., 1., 0., HBTK::Constants::pi()/2);
		}
		term_12 = HBTK::static_integrate(integrand, points, weights, n_pts);
		term_1  = term_11 * term_12;

		term_21 = -1. * HBTK::Constants::i() * (omega / U);
		term_22 = heave_added_mass / wing.area();
		term_2  = term_21 * term_22;

		assert(HBTK::check_finite(term_1));
		assert(HBTK::check_finite(term_2));
		return term_1 + term_2;
	}


	std::complex<double> Sclavounos1987::compute_lift_coeff_j5() 
	{
		assert(omega != 0.);
		assert(j == 5);
		std::complex<double> term_1, term_11, term_12;

		auto integrand = [&](double theta) -> std::complex<double> {
			auto y = wing.semispan() * cos(theta);
			auto semichord = wing.semichord(y);
			auto F_5 = F(y);
			auto C = Common::theodorsen_function((omega / U) * semichord);
			auto circulation = semichord + (2 * U / (HBTK::Constants::i() * omega) + 2. * F_5);
			auto added_mass = semichord * (1. + HBTK::Constants::i() * omega * F_5 / U);
			return HBTK::Constants::pi() * semichord * sin(theta) * (C * circulation + added_mass);
		};

		term_11 = 4. * wing.semispan() / wing.area();
		const int n_pts = 40;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			HBTK::linear_remap(points[idx], weights[idx], -1., 1., 0., HBTK::Constants::pi() / 2);
		}
		term_12 = HBTK::static_integrate(integrand, points, weights, n_pts);
		term_1 = term_11 * term_12;

		assert(HBTK::check_finite(term_1));
		return term_1;
	}


	double  Sclavounos1987::elliptic_added_mass_coefficient() {
		return mFlow::elliptic_added_mass_coefficient(wing.span, wing.chord(0));
	}


	double Sclavounos1987::rectangular_added_mass_coefficient()
	{
		return mFlow::rectangular_added_mass_coefficient(wing.span, wing.chord(0));
	}


	std::complex<double> Sclavounos1987::P(double y)
	{
		double term_1, term_2;

		// On the real term we're using integration by parts.
		auto integrand_1 = [y](const double t) -> double{
			return -y * exp(-y*t) * (asin(1./t) + sqrt(t*t - 1.) - t);
		};
		auto integrand_2 = [y](const double t) -> double{
			return exp(-y*t) * (sqrt(1 - t*t) - 1) / t;
		};
		
		const int n_points = 20;
		std::array<double, n_points> points1, weights1, points2, weights2;
		HBTK::gauss_legendre(points1, weights1);
		points2 = points1;
		weights2 = weights1;
		for (int idx = 0; idx < n_points; idx++) {
			HBTK::exponential_remap(points1[idx], weights1[idx], 1.);
			HBTK::linear_remap(points2[idx], weights2[idx], -1., 1., 0., 1.);
		}

		term_1 = - exp(-y) * (asin(1.) - 1)
			- HBTK::static_integrate(integrand_1, points1, weights1, n_points);
		term_2 = HBTK::static_integrate(integrand_2, points2, weights2, n_points);

		assert(HBTK::check_finite(term_1));
		assert(HBTK::check_finite(term_2));
		return term_1 + HBTK::Constants::i() * term_2;
	}

	double elliptic_added_mass_coefficient(double span, double max_chord)
	{
		// Reference: http://brennen.caltech.edu/fluidbook/basicfluiddynamics/unsteadyflows/addedmass/valuesoftheaddedmass.pdf
		// (Unable to find any closed form solution (even in Hydrodynamics(Lamb) 2nd Ed.)
		// Integration of the flat plate term in eq4.4 results in the equivalent of AR=1, without the correction
		// for the elliptic nature of the plate.
		double a = span / 2;
		double b = max_chord / 2;
		if (a < b) { std::swap(a, b); }
		double added_mass = a * b * b * HBTK::Constants::pi() * 4. / 3.;
		double ratio = a / b;
		// Linear interpolation of correction for ellipse from circle.
		std::vector<double> known_ratios = { 1., 1.5, 2., 3., 4., 6., 8.19, 10.34, 14.30, 10000. }; // To inf really.
		std::vector<double> known_coeffs = { 0.637, 0.748, 0.826, 0.900, 0.933, 0.964, 0.978, 0.985, 0.991, 1.000 };

		added_mass *= 2 * HBTK::linear_interpolate(known_ratios, known_coeffs, ratio);
		return added_mass;
	}

	double rectangular_added_mass_coefficient(double span, double chord)
	{
		// Reference: http://brennen.caltech.edu/fluidbook/basicfluiddynamics/unsteadyflows/addedmass/valuesoftheaddedmass.pdf
		// (Unable to find any closed form solution (even in Hydrodynamics(Lamb) 2nd Ed.)
		double a = span;
		double b = chord;
		if (a < b) { std::swap(a, b); }
		double added_mass = a * b * b * HBTK::Constants::pi() / 4.;
		double ratio = a / b;
		// Linear interpolation of correction for ellipse from circle.
		std::vector<double> known_ratios = { 1., 1.5, 2., 3., 10000. }; // To inf really.
		std::vector<double> known_coeffs = { 0.478, 0.680, 0.840, 1.000, 1.000 };
		added_mass *= 2 * HBTK::linear_interpolate(known_ratios, known_coeffs, ratio);
		return added_mass;
	}

}
