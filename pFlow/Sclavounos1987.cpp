#include "stdafx.h"
#include "Sclavounos1987.h"
/*////////////////////////////////////////////////////////////////////////////
Guermond1991.cpp

Equations from the paper "An unsteady lifting-line theory", P.D. Sclavounos,
1987, Journal of Engineering Mathmatics

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

#include <Eigen/LU>

#include "Common.h"

namespace mFlow {
	Sclavounos1987::Sclavounos1987()
	{
	}


	Sclavounos1987::~Sclavounos1987()
	{
	}


	void Sclavounos1987::compute_collocation_points()
	{
		// After Eq 7.1 states these as local maxima of sin((N+1)theta). This doesn't work.
		m_collocation_points.resize(number_of_terms);
		const auto hpi = HBTK::Constants::pi() / 2;
		for (auto idx = 0; idx < number_of_terms; idx++) {
			m_collocation_points[idx] = (hpi + idx * HBTK::Constants::pi()) / (2 * number_of_terms);
			//m_collocation_points[idx] = hpi / (idx + 1);
			//m_collocation_points[idx] = hpi / (2*idx + 1);
			assert(m_collocation_points[idx] <= HBTK::Constants::pi() && m_collocation_points[idx] >= 0.);
		}
		return;
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
		auto E_1_term = -HBTK::Constants::i() * nu * Common::exponential_int_E1(nu * abs(y));

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
		auto quad = get_split_quad(num_quad_points, theta_sing);

		// We're using the singularity subtraction method.
		auto singularity_coefficient = dsintheta_dtheta(theta_sing, k) * K_term1_numerator(0);

		auto numerical_integrand = [&](double theta0){
			auto eta = wing.semispan() * cos(theta0);
			auto singular_part = K_term1_singularity(y_position - eta);
			auto nonsingular_K = K_term1_numerator(y_position - eta);
			auto gamma_dtheta = dsintheta_dtheta(theta0, k);

			auto singularity_subtraction = (nonsingular_K * gamma_dtheta - singularity_coefficient);
			return singular_part * singularity_subtraction;
		};

		auto integral = HBTK::static_integrate(numerical_integrand, std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size())
			+ HBTK::static_integrate(numerical_integrand, std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size())
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

		auto integral_coefficient = HBTK::Constants::i() * omega / (2 * U);

		// The non-singular part of the integrand
		auto non_singular = [&](double theta)->double {
			return (2 * k + 1) * cos((2 * k + 1)*theta) / (wing.semispan() * sin(theta));
		};
		// The singular part of the integrand
		auto singular = [&](double theta)->std::complex<double> {
			return wing.semispan() * sin(theta) * (theta > theta_sing ? 1. : -1.)
				* Common::exponential_int_E1(wing.semispan() * abs(cos(theta_sing) - cos(theta)));
		};
		// The singular part of the integrand evaluated.
		auto singular_integral = wing.semispan() * ((cos(theta_sing)+1.) * Common::exponential_int_E1(wing.semispan() * abs(cos(theta_sing)+1.))
			+ (cos(theta_sing) -1.) * Common::exponential_int_E1(wing.semispan() * abs(cos(theta_sing)-1.))) 
			+ exp(wing.semispan()*(cos(theta_sing) - 1.)) - exp(-wing.semispan()*(cos(theta_sing) + 1.));

		// The singularity subtraction method is used here.
		auto ssm_variable = non_singular(theta_sing);
		auto numerical_integrand = [&](double theta) {
			auto singular_var = singular(theta);
			auto non_singular_var = non_singular(theta);
			auto singularity_subtraction = non_singular_var - ssm_variable;
			auto integrand = singular_var * singularity_subtraction;
			return integrand;
		};

		auto quad = get_split_quad(100, theta_sing);
		auto int_lower = HBTK::static_integrate(numerical_integrand, std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size());
		auto int_upper = HBTK::static_integrate(numerical_integrand, std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size());
		auto singular_interal = ssm_variable * singular_integral;

		auto complete_integral = integral_coefficient * (int_lower + int_upper + singular_integral);
		assert(HBTK::check_finite(complete_integral));

		return complete_integral;

		// We want to split our integral to remove the abs(y-eta) and sign(y-eta) parts.
		// Additionally, we'll using integration by parts followed by the singularity
		// subtraction method. All integrating with respect to theta.
		/*
		// I = int(u dv) = uv - int(du v).
		// dv = (2k+1) * cos((2k+1)theta), u = - 0.5 * sgn(y-eta) * v * i * E_1(v * abs(y-eta)) 

		// integral of the fourier part to v
		auto fourier_term = [&](double theta_0) {
			return sin((2 * k + 1) * theta_0);
		};

		// Derivative of the K terms (ie. du)
		auto K_term2_dtheta_singular = [&](double ypos, double eta) {
			return  wing.semispan() / (ypos - eta);
		};
		auto K_term2_dtheta_singular_integral = [&](double theta) {
			return 0.0; // From 0 to pi this is a Glauert Integral.
		};
		auto K_term2_dtheta_numerator = [&](double ypos, double eta) {
			return HBTK::Constants::i() * (omega / U) * 0.5 * sin(acos(eta / wing.semispan())) 
				* exp(-(omega / U) * abs(ypos - eta));
		};

		// compute uv:
		auto uv = HBTK::Constants::i() * (omega / U) * Common::Exponential_int_E1(0.) * sin((2 * k + 1)*theta_sing);
		
		// compute int u dv		
		auto ssm_coeff = fourier_term(theta_sing) * K_term2_dtheta_numerator(y, y);
		auto integrand = [&](double theta_0) {
			auto eta = wing.semispan() * cos(theta_0);
			auto singular_bit = K_term2_dtheta_singular(y, eta);
			auto nonsingular = fourier_term(theta_0) * K_term2_dtheta_numerator(y, eta);
			return singular_bit * (nonsingular - ssm_coeff);
		};
		auto quad = get_split_quad(40, theta_sing);
		auto intudv1 = HBTK::static_integrate(integrand, std::get<0>(quad), std::get<1>(quad), std::get<0>(quad).size());
		auto intudv2 = HBTK::static_integrate(integrand, std::get<2>(quad), std::get<3>(quad), std::get<2>(quad).size());
		auto intudv3 = ssm_coeff * K_term2_dtheta_singular_integral(theta_sing); // The integral of the singular bit is 0.
		auto intudv = intudv1 + intudv2 + intudv3;

		assert(HBTK::check_finite(intudv));
		return -(uv + intudv);
		*/
	}


	std::complex<double> Sclavounos1987::integrate_gammaprime_K_term3(double y, int k)
	{
		// Integrating the P term of of int ( Gamma'(eta) K(y-eta) deta).
		// We're integrating in [0, singularity], [singularity, pi]
		auto theta_sing = acos(y / wing.semispan());
		auto quad = get_split_quad(40, theta_sing); // lower_points, lower_weights, upper_points, upper_weights
		auto integrand = [&](double theta0) {
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
		std::complex<double> F_res = 0.;
		for (int idx = 0; idx < (int)m_solution.size(); idx++) {
			F_res += integrate_gammaprime_K(y, idx) * m_solution[idx];
		}
		F_res = F_res * -1. / (HBTK::Constants::i() * omega * 2. * HBTK::Constants::pi());
		assert(HBTK::check_finite(F_res));
		return F_res;
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
		Sclavounos1987::get_split_quad(int total_pts, double split_theta)
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

		compute_collocation_points();

		// Allocation of matrices.
		m_solution.resize(number_of_terms);
		Eigen::MatrixXd gamma_matrix; // The matrix that directly gives vorticities at given points.
		// The matrix that representsthe integro-differential part.
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> integ_diff_matrix;
		// For the probem in for Ax = b, solve for x, this is A.
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> LHS_matrix;
		// RHS vector
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> RHS_vector;
		gamma_matrix.resize(number_of_terms, number_of_terms);
		integ_diff_matrix.resize(number_of_terms, number_of_terms);
		m_gammaprime_K_matrix = integ_diff_matrix;
		LHS_matrix.resize(number_of_terms, number_of_terms);
		RHS_vector.resize(number_of_terms);

		// Compute gamma_matrix;
		for (int i = 0; i < number_of_terms; i++) {
			for (int j = 0; j < number_of_terms; j++) {
				gamma_matrix(i, j) = sin((2 * j + 1) * m_collocation_points[i]);
			}
		}

		// Compute integ_diff_matrix
		for (int i = 0; i < number_of_terms; i++) {

			auto y_position = wing.semispan() * cos(m_collocation_points[i]);
			std::complex<double> ext_coeff;
			if (omega != 0) { 
				ext_coeff = d_3(y_position)	/ (2 * HBTK::Constants::pi() * omega * HBTK::Constants::i());
			} else {
				ext_coeff = - wing.semichord(y_position);
			}

			for (int j = 0; j < number_of_terms; j++) 
			{
				auto integral = integrate_gammaprime_K(y_position, j);
				m_gammaprime_K_matrix(i, j) = integral;
				integ_diff_matrix(i, j) = ext_coeff * integral;
			} // End For in j
		} // End for in i

		LHS_matrix = gamma_matrix - integ_diff_matrix;
		
		// Generate RHS vector:
		if (j == 3) {
			for (int i = 0; i < number_of_terms; i++) {
				auto y = wing.semispan()*cos(m_collocation_points[i]);
				RHS_vector[i] = d_3(y);
			}
		}
		else if (j == 5) {
			for (int i = 0; i < number_of_terms; i++) {
				auto y = wing.semispan()*cos(m_collocation_points[i]);
				RHS_vector[i] = d_5(y) - U * d_3(y) / (HBTK::Constants::i() * omega);
			}
		}

		assert(abs(LHS_matrix.determinant()) >= 1e-5);  
		m_solution = LHS_matrix.lu().solve(RHS_vector); // Solve matrix problem using LU decomposition.
		return;
	}


	std::complex<double> Sclavounos1987::get_solution_vorticity(double y)
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


	std::complex<double> Sclavounos1987::compute_lift_coeff(double heave_added_mass)
	{
		assert(omega != 0.);
		std::complex<double> term_1, term_2,
			term_11, term_12, term_21, term_22;

		auto integrand = [&](double theta) -> std::complex<double> {
			auto y = wing.semispan() * cos(theta);
			auto semichord = wing.semichord(y);
			auto F_3 = F(y);
			auto C = Common::theodorsen_function((omega / U) * semichord);
			return C * semichord * F_3 * sin(theta);
		};

		term_11 = -4. / wing.area();
		const int n_pts = 60;
		std::array<double, n_pts> points, weights;
		HBTK::gauss_legendre(points, weights);
		for (int idx = 0; idx < n_pts; idx++) {
			HBTK::linear_remap(points[idx], weights[idx], -1., 1., 0., HBTK::Constants::pi()/2);
		}
		term_12 = wing.semispan() * 2.0 * HBTK::static_integrate(integrand, points, weights, n_pts);
		term_1 = term_11 * term_12;

		term_21 = HBTK::Constants::i() * (omega / U);
		term_22 = heave_added_mass / wing.area();
		term_2 = term_21 * term_22;

		assert(HBTK::check_finite(term_1));
		assert(HBTK::check_finite(term_2));
		return term_1 - term_2;
	}


	double  Sclavounos1987::get_elliptic_added_mass_coefficient() {
		// Reference: http://brennen.caltech.edu/fluidbook/basicfluiddynamics/unsteadyflows/addedmass/valuesoftheaddedmass.pdf
		// (Unable to find any closed form solution (even in Hydrodynamics(Lamb) 2nd Ed.)
		// Integration of the flat plate term in eq4.4 results in the equivalent of AR=1, without the correction
		// for the elliptic nature of the plate.
		double a = wing.semispan();
		double b = wing.semichord(0.);
		if (a < b) { std::swap(a, b); }
		double added_mass = a * b * b * HBTK::Constants::pi() * 4. / 3.;
		double ratio = a / b;
		// Linear interpolation of correction for ellipse from circle.
		std::vector<double> known_ratios = { 1., 1.5, 2., 3., 4., 6., 8.19, 10.34, 14.30, 10000. }; // To inf really.
		std::vector<double> known_coeffs = { 0.637, 0.748, 0.826, 0.900, 0.933, 0.964, 0.978, 0.985, 0.991, 1.000 };
		double coeff = 0;
		int lower_known = -1;
		for (int idx = 0; idx < (int)known_ratios.size(); idx++) {
			if (ratio >= known_ratios[idx]) {
				lower_known = idx;
			}
		}
		if (lower_known == 0) { coeff = known_ratios[0]; } // == 1.

		double fraction = (ratio - known_ratios[lower_known]) / (known_ratios[lower_known + 1] - known_ratios[lower_known]);
		coeff = (known_coeffs[lower_known + 1] - known_coeffs[lower_known]) * fraction + known_coeffs[lower_known];
		added_mass *= 2 * coeff;
		return added_mass;
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

}
