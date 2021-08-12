#include "libphysica/Integration.hpp"

#include <algorithm>
#include <iostream>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

namespace libphysica
{
using namespace boost::math::quadrature;

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
void Check_Integration_Limits(double& a, double& b, double& sign)
{
	if(a > b)
	{
		std::cerr << "Warning in libphysica::Integrate(): From the integral from a to b, a>b (a = " << a << ", b = " << b << "). Sign will get swapped." << std::endl;
		std::swap(a, b);
		sign = -1.0;
	}
}

double Integrate(std::function<double(double)> func, double a, double b, double epsilon, int maxRecursionDepth)
{
	double sign = +1.0;
	if(a == b)
		return 0.0;
	else
		Check_Integration_Limits(a, b, sign);
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

// 1.2 1D integration with boost functions
double Integrate(std::function<double(double)> func, double a, double b, const std::string& method)
{
	double sign = 1.0;
	if(a == b)
		return 0.0;
	else
		Check_Integration_Limits(a, b, sign);
	if(method == "Trapezoidal")
		return sign * trapezoidal(func, a, b);
	else if(method == "Gauss-Legendre")
		return sign * gauss<double, 30>::integrate(func, a, b);
	else if(method == "Gauss-Kronrod")
		return sign * gauss_kronrod<double, 31>::integrate(func, a, b, 5, 1e-9);
	else
	{
		std::cerr << "Error in libphysica::Integrate(): Method " << method << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// 2. Multidimensional integration
// 2.1 Multidimensional integration via nesting 1D integration
double Integrate_2D(std::function<double(double, double)> func, double x1, double x2, double y1, double y2, const std::string& method)
{
	auto integrand_x = [&func, y1, y2, method](double x) {
		auto integrand_y = [&func, x](double y) {
			return func(x, y);
		};
		return Integrate(integrand_y, y1, y2, method);
	};
	return Integrate(integrand_x, x1, x2, method);
}

double Integrate_3D(std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2, const std::string& method)
{
	auto integrand_x = [&func, y1, y2, z1, z2, method](double x) {
		auto integrand_y = [&func, z1, z2, x, method](double y) {
			auto integrand_z = [&func, x, y, method](double z) {
				return func(x, y, z);
			};
			return Integrate(integrand_z, z1, z2, method);
		};
		return Integrate(integrand_y, y1, y2, method);
	};
	return Integrate(integrand_x, x1, x2, method);
}

double Integrate_3D(std::function<double(Vector)> func, double r1, double r2, double costheta_1, double costheta_2, double phi_1, double phi_2, const std::string& method)
{
	auto integrand = [&func](double r, double cos_theta, double phi) {
		Vector rVec = Spherical_Coordinates(r, acos(cos_theta), phi);
		return r * r * func(rVec);
	};
	return Integrate_3D(integrand, r1, r2, costheta_1, costheta_2, phi_1, phi_2, method);
}

// 2.2 Monte Carlo Integration
// Reference: Some of these functions are taken from http://numerical.recipes/webnotes/nr3web9.pdf

// Utility routine used by vegas, to rebin a vector of densities contained in row j of xi into new bins defined by a vector r.
void Rebin(const double rc, const int nd, std::vector<double>& r, std::vector<double>& xin, libphysica::Matrix& xi, const int j)
{
	int i, k = 0;
	double dr = 0.0, xn = 0.0, xo = 0.0;

	for(i = 0; i < nd - 1; i++)
	{
		while(rc > dr)
			dr += r[(++k) - 1];
		if(k > 1)
			xo = xi[j][k - 2];
		xn = xi[j][k - 1];
		dr -= rc;
		xin[i] = xn - (xn - xo) * dr / r[k - 1];
	}
	for(i = 0; i < nd - 1; i++)
		xi[j][i] = xin[i];
	xi[j][nd - 1] = 1.0;
}

double Integrate_MC_Vegas(std::function<double(std::vector<double>&, const double)> func, std::vector<double>& region, const int init, const int ncall, const int itmx, const int nprn)
{
	double integral, chi2a, standard_deviation;
	// Best make everything static, allowing restarts.
	static const int NDMX = 50, MXDIM = 10;
	static const double ALPH = 1.5, TINY = 1.0e-30;
	static int i, it, j, k, mds, nd, ndo, ng, npg;
	static double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti;
	static double tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt;
	static std::vector<int> ia(MXDIM), kg(MXDIM);
	static std::vector<double> dt(MXDIM), dx(MXDIM), r(NDMX), x(MXDIM), xin(NDMX);
	static libphysica::Matrix d(NDMX, MXDIM), di(NDMX, MXDIM), xi(MXDIM, NDMX);

	// Initialize  captive, static random number generator
	std::random_device rd;
	std::mt19937 PRNG(rd());

	int ndim = region.size() / 2;
	if(init <= 0)
	{
		mds = ndo = 1;
		for(j = 0; j < ndim; j++)
			xi[j][0] = 1.0;
	}
	if(init <= 1)
		si = swgt = schi = 0.0;
	if(init <= 2)
	{
		nd = NDMX;
		ng = 1;
		if(mds != 0)
		{
			ng	= int(pow(ncall / 2.0 + 0.25, 1.0 / ndim));
			mds = 1;
			if((2 * ng - NDMX) >= 0)
			{
				mds = -1;
				npg = ng / NDMX + 1;
				nd	= ng / npg;
				ng	= npg * nd;
			}
		}
		for(k = 1, i = 0; i < ndim; i++)
			k *= ng;
		npg	  = std::max(int(ncall / k), 2);
		calls = double(npg) * double(k);
		dxg	  = 1.0 / ng;
		for(dv2g = 1, i = 0; i < ndim; i++)
			dv2g *= dxg;
		dv2g = calls * dv2g * calls * dv2g / npg / npg / (npg - 1.0);
		xnd	 = nd;
		dxg *= xnd;
		xjac = 1.0 / calls;
		for(j = 0; j < ndim; j++)
		{
			dx[j] = region[j + ndim] - region[j];
			xjac *= dx[j];
		}
		if(nd != ndo)
		{
			for(i = 0; i < std::max(nd, ndo); i++)
				r[i] = 1.0;
			for(j = 0; j < ndim; j++)
				Rebin(ndo / xnd, nd, r, xin, xi, j);
			ndo = nd;
		}
		if(nprn >= 0)
		{
			std::cout << " Input parameters for vegas";
			std::cout << "  ndim= " << std::setw(4) << ndim;
			std::cout << "  ncall= " << std::setw(8) << calls << std::endl;
			std::cout << std::setw(34) << "  it=" << std::setw(5) << it;
			std::cout << "  itmx=" << std::setw(5) << itmx << std::endl;
			std::cout << std::setw(34) << "  nprn=" << std::setw(5) << nprn;
			std::cout << "  ALPH=" << std::setw(9) << ALPH << std::endl;
			std::cout << std::setw(34) << "  mds=" << std::setw(5) << mds;
			std::cout << "  nd=" << std::setw(5) << nd << std::endl;
			for(j = 0; j < ndim; j++)
			{
				std::cout << std::setw(30) << " x1[" << std::setw(2) << j;
				std::cout << "]= " << std::setw(11) << region[j] << " xu[";
				std::cout << std::setw(2) << j << "]= ";
				std::cout << std::setw(11) << region[j + ndim] << std::endl;
			}
		}
	}
	for(it = 0; it < itmx; it++)
	{
		ti = tsi = 0.0;
		for(j = 0; j < ndim; j++)
		{
			kg[j] = 1;
			for(i = 0; i < nd; i++)
				d[i][j] = di[i][j] = 0.0;
		}
		for(;;)
		{
			fb = f2b = 0.0;
			for(k = 0; k < npg; k++)
			{
				wgt = xjac;
				for(j = 0; j < ndim; j++)
				{
					xn	  = (kg[j] - libphysica::Sample_Uniform(PRNG)) * dxg + 1.0;
					ia[j] = std::max(std::min(int(xn), NDMX), 1);
					if(ia[j] > 1)
					{
						xo = xi[j][ia[j] - 1] - xi[j][ia[j] - 2];
						rc = xi[j][ia[j] - 2] + (xn - ia[j]) * xo;
					}
					else
					{
						xo = xi[j][ia[j] - 1];
						rc = (xn - ia[j]) * xo;
					}
					x[j] = region[j] + rc * dx[j];
					wgt *= xo * xnd;
				}
				f  = wgt * func(x, wgt);
				f2 = f * f;
				fb += f;
				f2b += f2;
				for(j = 0; j < ndim; j++)
				{
					di[ia[j] - 1][j] += f;
					if(mds >= 0)
						d[ia[j] - 1][j] += f2;
				}
			}
			f2b = sqrt(f2b * npg);
			f2b = (f2b - fb) * (f2b + fb);
			if(f2b <= 0.0)
				f2b = TINY;
			ti += fb;
			tsi += f2b;
			if(mds < 0)
			{
				for(j = 0; j < ndim; j++)
					d[ia[j] - 1][j] += f2b;
			}
			for(k = ndim - 1; k >= 0; k--)
			{
				kg[k] %= ng;
				if(++kg[k] != 1)
					break;
			}
			if(k < 0)
				break;
		}
		tsi *= dv2g;
		wgt = 1.0 / tsi;
		si += wgt * ti;
		schi += wgt * ti * ti;
		swgt += wgt;
		integral = si / swgt;
		chi2a	 = (schi - si * integral) / (it + 0.0001);
		if(chi2a < 0.0)
			chi2a = 0.0;
		standard_deviation = sqrt(1.0 / swgt);
		tsi				   = sqrt(tsi);
		if(nprn >= 0)
		{
			std::cout << " iteration no. " << std::setw(3) << (it + 1);
			std::cout << " : integral = " << std::setw(14) << ti;
			std::cout << " +/- " << std::setw(9) << tsi << std::endl;
			std::cout << " all iterations:  "
					  << " integral =";
			std::cout << std::setw(14) << integral << "+-" << std::setw(9) << standard_deviation;
			std::cout << " chi**2/IT n =" << std::setw(9) << chi2a << std::endl;
			if(nprn != 0)
			{
				for(j = 0; j < ndim; j++)
				{
					std::cout << " DATA FOR axis  " << std::setw(2) << j << std::endl;
					std::cout << "     X      delta i          X      delta i";
					std::cout << "          X       deltai" << std::endl;
					for(i = nprn / 2; i < nd - 2; i += nprn + 2)
					{
						std::cout << std::setw(8) << xi[j][i] << std::setw(12) << di[i][j];
						std::cout << std::setw(12) << xi[j][i + 1] << std::setw(12) << di[i + 1][j];
						std::cout << std::setw(12) << xi[j][i + 2] << std::setw(12) << di[i + 2][j];
						std::cout << std::endl;
					}
				}
			}
		}
		for(j = 0; j < ndim; j++)
		{
			xo		= d[0][j];
			xn		= d[1][j];
			d[0][j] = (xo + xn) / 2.0;
			dt[j]	= d[0][j];
			for(i = 2; i < nd; i++)
			{
				rc			= xo + xn;
				xo			= xn;
				xn			= d[i][j];
				d[i - 1][j] = (rc + xn) / 3.0;
				dt[j] += d[i - 1][j];
			}
			d[nd - 1][j] = (xo + xn) / 2.0;
			dt[j] += d[nd - 1][j];
		}
		for(j = 0; j < ndim; j++)
		{
			rc = 0.0;
			for(i = 0; i < nd; i++)
			{
				if(d[i][j] < TINY)
					d[i][j] = TINY;
				r[i] = pow((1.0 - d[i][j] / dt[j]) /
							   (log(dt[j]) - log(d[i][j])),
						   ALPH);
				rc += r[i];
			}
			Rebin(rc / xnd, nd, r, xin, xi, j);
		}
	}
	return integral;
}

double Integrate_MC_Brute_Force(std::function<double(std::vector<double>&, const double)> func, std::vector<double>& region, const int ncall)
{
	int dim = region.size() / 2.0;
	std::random_device rd;
	std::mt19937 PRNG(rd());

	double volume = 1.0;
	for(int i = 0; i < dim; i++)
		volume *= (region[i + dim] - region[i]);

	double sum = 0.0;
	// double sum_2 = 0.0;
	for(int i = 0; i < ncall; i++)
	{
		std::vector<double> args(dim);
		for(int j = 0; j < dim; j++)
			args[j] = region[j] + libphysica::Sample_Uniform(PRNG) * (region[j + dim] - region[j]);
		double fct = func(args, 0.0);
		sum += volume * fct;
		// sum_2 += volume * volume * fct * fct;
	}
	double integral = sum / ncall;
	// double standard_deviation = sqrt((sum_2 / ncall - integral * integral) / (ncall - 1.0));
	// double chi2a			  = 0.0;
	return integral;
}

double Integrate_MC(std::function<double(std::vector<double>&, const double)> func, std::vector<double>& region, const int ncalls, const std::string& method)
{
	if(method == "Brute force")
		return Integrate_MC_Brute_Force(func, region, ncalls);
	else if(method == "Vegas")
	{
		int init	   = 0;
		const int itmx = 5;
		const int nprn = -1;
		return Integrate_MC_Vegas(func, region, init, ncalls, itmx, nprn);
	}
	else
	{
		std::cerr << "Error in libphysica::Integrate_MC(): Method " << method << " not recognized." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace libphysica