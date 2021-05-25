#include "libphysica/Integration.hpp"

#include <cmath>
#include <iostream>

#include "libphysica/Special_Functions.hpp"

namespace libphysica
{

// 1. One-dimensional MC integration
// 1.1 One-dimensional integration via adaptive Simpson method

//Function to return a reasonable precision.
double Find_Epsilon(std::function<double(double)> func, double a, double b, double precision)
{
	double c	   = (a + b) / 2;
	double h	   = b - a;
	double fa	   = func(a);
	double fb	   = func(b);
	double fc	   = func(c);
	double S	   = (h / 6) * (fa + 4 * fc + fb);
	double epsilon = precision * S;
	return epsilon;
}

double Adaptive_Simpson_Integration(std::function<double(double)> func, double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom, bool& warning)
{
	double c	  = (a + b) / 2;
	double h	  = b - a;
	double d	  = (a + c) / 2;
	double e	  = (b + c) / 2;
	double fd	  = func(d);
	double fe	  = func(e);
	double Sleft  = (h / 12) * (fa + 4 * fd + fc);
	double Sright = (h / 12) * (fc + 4 * fe + fb);
	double S2	  = Sleft + Sright;
	if(bottom <= 0 || fabs(S2 - S) <= 15 * epsilon)	  //15 due to error analysis
	{
		if(bottom <= 0 && fabs(S2 - S) > 15 * epsilon)
			warning = true;
		return S2 + (S2 - S) / 15;
	}
	else
	{
		return Adaptive_Simpson_Integration(func, a, c, epsilon / 2, Sleft, fa, fc, fd, bottom - 1, warning) + Adaptive_Simpson_Integration(func, c, b, epsilon / 2, Sright, fc, fb, fe, bottom - 1, warning);
	}
}

//Recursive functions for one-dimensional integration
double Integrate(std::function<double(double)> func, double a, double b, double epsilon, int maxRecursionDepth)
{
	int sign = +1;
	if(a == b)
		return 0.0;
	else if(a > b)
	{
		double aux = a;
		a		   = b;
		b		   = aux;
		sign	   = -1;
	}
	double c	  = (a + b) / 2;
	double h	  = b - a;
	double fa	  = func(a);
	double fb	  = func(b);
	double fc	  = func(c);
	double S	  = (h / 6) * (fa + 4 * fc + fb);
	bool warning  = false;
	double result = Adaptive_Simpson_Integration(func, a, b, fabs(epsilon), S, fa, fb, fc, maxRecursionDepth, warning);
	if(warning)
	{
		std::cout << "Warning in libphysica::Integrate(): Numerical integration on the interval (" << a << "," << b << ") did not converge to the desired precision." << std::endl;
		std::cout << "\tDesired precision: " << Round(fabs(epsilon)) << " Result: " << Round(result) << std::endl;
	}
	if(std::isnan(result))
		std::cout << "Warning in libphysica::Integrate(): Result is nan." << std::endl;
	else if(std::isinf(result))
		std::cout << "Warning in libphysica::Integrate(): Result is inf." << std::endl;
	return sign * result;
}
// 2. Multidimensional MC integration
void Integrate_MC_Vegas(std::function<double(std::vector<double>&)> fxn, std::vector<double>& region, const int init, const int ncall, const int itmx, const int nprn, double& sd, double& chi2a)
{
}

}	// namespace libphysica