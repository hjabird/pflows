#include "CanonicalFuntions.h"

#include <cassert>
#include <cmath>

#include <HBTK/Checks.h>

mFlow::Harmonic::Harmonic(double omega, double amplitude, double phase_offset)
	: amplitude(amplitude),
	angular_vel(omega),
	phase_diff(phase_offset)
{
}

double mFlow::Harmonic::f(double x)
{
	return amplitude * sin(angular_vel * x + phase_diff);
}

double mFlow::Harmonic::dfdx(double x)
{
	return angular_vel * amplitude * cos(angular_vel * x + phase_diff);
}

mFlow::EldredgeSmoothRamp::EldredgeSmoothRamp(double t_rampup_start, double t_rampup_end,
	double t_rampdown_start, double t_rampdown_end,
	double smoothing, double amplitude, double chord, double free_stream_vel)
	: t1(t_rampup_start),
	t2(t_rampup_end),
	t3(t_rampdown_start),
	t4(t_rampdown_end),
	a(smoothing),
	A(amplitude),
	c(chord),
	U(free_stream_vel)
{
	assert(t1 < t2);
	assert(t2 <= t3);
	assert(t3 < t4);
	assert(a > 0);
	assert(chord > 0);
	assert(U > 0);
	assert(HBTK::check_finite(A));
	assert(HBTK::check_finite(c));
	assert(HBTK::check_finite(U));
}

double mFlow::EldredgeSmoothRamp::f(double t)
{
	double valG = G(t);
	// I need a numerical method to do the following properly:
	double maxG = G((t3 + t2) / 2); 
	return A * valG / maxG;
}

double mFlow::EldredgeSmoothRamp::dfdx(double t)
{
	double valdG = dGdt(t);
	// I need a numerical method to do the following properly:
	double maxG = G((t3 + t2) / 2);
	return A * valdG / maxG;
}

inline double mFlow::EldredgeSmoothRamp::G(double t)
{
	return log(eld_cosh(t, t1) * eld_cosh(t, t4) / (eld_cosh(t, t2) * eld_cosh(t, t3)));
}

inline double mFlow::EldredgeSmoothRamp::dGdt(double t)
{
	double term_1, term_2, term_2n, term_2d,
		term_2n1, term_2n2, term_2n3, term_2n4;
	term_2n1 = eld_cosh(t, t1) * eld_sinh(t, t4) /
		(eld_cosh(t, t2) * eld_cosh(t, t3));
	term_2n2 = eld_cosh(t, t1) * eld_sinh(t, t3) * eld_cosh(t, t4) /
		(eld_cosh(t, t2) * pow(eld_cosh(t, t3), 2));
	term_2n3 = eld_cosh(t, t1) * eld_sinh(t, t2) * eld_cosh(t, t4) /
		(pow(eld_cosh(t, t2), 2) * eld_cosh(t, t3));
	term_2n4 = eld_sinh(t, t1) * eld_cosh(t, t4) /
		(eld_cosh(t, t2) * eld_cosh(t, t3));
	term_2n = term_2n1 - term_2n2 - term_2n3 + term_2n4;
	term_2d = eld_cosh(t, t1) * eld_cosh(t, t4);
	term_2 = term_2n / term_2d;
	term_1 = eld_cosh(t, t2) * eld_cosh(t, t3) * U * a / c;
	return term_1 * term_2;
}

inline double mFlow::EldredgeSmoothRamp::eld_cosh(double t, double t_ref)
{
	return cosh((U * a * (t - t_ref)) / c);
}

inline double mFlow::EldredgeSmoothRamp::eld_sinh(double t, double t_ref)
{
	return sinh((U * a * (t - t_ref)) / c);
}

double mFlow::CanonicalFunction::f(double x)
{
	return 0.0;
}

double mFlow::CanonicalFunction::dfdx(double x)
{
	return 0.0;
}

mFlow::ConstantValue::ConstantValue(double value)
	: value(value)
{
	assert(HBTK::check_finite(value));
}

double mFlow::ConstantValue::f(double x)
{
	return value;
}

double mFlow::ConstantValue::dfdx(double x)
{
	return 0.0;
}

mFlow::CanonicalFunctionSum::CanonicalFunctionSum(std::unique_ptr<CanonicalFunction> summand_1,
	std::unique_ptr<CanonicalFunction> summand_2)
	: summand_1(std::move(summand_1)),
	summand_2(std::move(summand_2))
{
}

	double mFlow::CanonicalFunctionSum::f(double x)
{
	return summand_1->f(x) + summand_2->f(x);
}

double mFlow::CanonicalFunctionSum::dfdx(double x)
{
	return summand_1->dfdx(x) + summand_2->dfdx(x);
}
