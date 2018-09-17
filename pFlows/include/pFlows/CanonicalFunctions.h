#pragma once
/*////////////////////////////////////////////////////////////////////////////
CanonicalFunctions.h

Functions that are often encountered where one does not wish to type out the
derivative etc etc every time one encounters them.

Copyright 2018 HJA Bird

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

#include <cassert>
#include <memory>

namespace mFlow {
	// Base class. Returns a zero for anything.
	class CanonicalFunction {
	public:
		CanonicalFunction() = default;
		virtual ~CanonicalFunction() = default;
		virtual double f(double x);
		virtual double dfdx(double x);
	};

	// Harmonic oscillation - amplitude * sin( omega * t + phase_offset)
	class Harmonic : 
		public CanonicalFunction
	{
	public:
		Harmonic(double omega, double amplitude, double phase_offset);

		virtual double f(double x) override;
		virtual double dfdx(double x) override;

		double amplitude;	// The difference between mean and peak
		double angular_vel;	// sin( OMEGA * t) where OMEGA is angular_vel
		double phase_diff;	// sin(w*t + PHI) where PHI is phase difference
	};


	// A ramp hold return function defined by Eldredge
	// f(t) = Amplitude * G(t) / max(G(t)) where G(t) is Eldredge's smooth ramp
	// See Ramesh 2013 for example.
	class EldredgeSmoothRamp : 
		public CanonicalFunction
	{
	public:
		EldredgeSmoothRamp(double t_rampup_start, double t_rampup_end, 
			double t_rampdown_start, double t_rampdown_end, 
			double smoothing, double amplitude, double chord,
			double free_stream_vel);

		virtual double f(double t) override;
		virtual double dfdx(double t) override;

		double a;	// Defines smoothing at corners
		double t1;	// Time from 0 until start of ramp motion
		double t2;	// t1 + duration of upstroke (end of upstroke)
		double t3;	// t2 + hold time at max
		double t4;	// t3 + duration of downstroke
		double A;	// Amplitude
		double c;	// Chord
		double U;	// Free stream vel

	private:
		inline double G(double t);
		inline double dGdt(double t);
		inline double eld_cosh(double t, double t_ref);
		inline double eld_sinh(double t, double t_ref);

	};

	// A static value.
	class ConstantValue : 
		public CanonicalFunction
	{
	public:
		ConstantValue(double value);
		virtual double f(double x) override;
		virtual double dfdx(double x) override;
		double value;
	};

	// Sum of two functions.
	class CanonicalFunctionSum : 
		public CanonicalFunction
	{
	public:
		CanonicalFunctionSum(std::unique_ptr<CanonicalFunction> summand_1, 
			std::unique_ptr<CanonicalFunction> summand_2);
		virtual double f(double x) override;
		virtual double dfdx(double x) override;
		std::unique_ptr<CanonicalFunction> summand_1, summand_2;
	};
}