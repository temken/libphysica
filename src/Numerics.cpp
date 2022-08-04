// Disclaimer:
// Some of the function implementations were made with the help of the
//"Numerical Recipes 3rd Edition: The Art of Scientific Computing"
// by William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery

#include "libphysica/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

namespace libphysica
{

// 1. Interpolation
// 1.1 One-dimensional interpolation
void Interpolation::Compute_Steffen_Coefficients()
{
	a.clear();
	b.clear();
	c.clear();
	d.clear();

	// Compute the Steffen coefficients for the interpolation
	// 1. h and s.
	std::vector<double> h(N - 1), s(N - 1);
	for(unsigned int i = 0; i < N - 1; i++)
	{
		h[i] = x_values[i + 1] - x_values[i];
		s[i] = (function_values[i + 1] - function_values[i]) / h[i];
	}

	// 2. p and dy
	std::vector<double> dy(N), p(N);
	for(unsigned int i = 0; i < N; i++)
	{
		// First point
		if(i == 0)
		{
			p[i]  = s[i] * (1.0 + h[i] / (h[i] + h[i + 1])) - s[i + 1] * h[i] / (h[i] + h[i + 1]);
			dy[i] = (Sign(p[i]) + Sign(s[i])) * std::min(1.0 * fabs(s[i]), 0.5 * fabs(p[i]));
		}
		// Last point
		else if(i == N - 1)
		{
			p[i]  = s[i - 1] * (1.0 + h[i - 1] / (h[i - 1] + h[i - 2])) - s[i - 2] * h[i - 1] / (h[i - 1] + h[i - 2]);
			dy[i] = (Sign(p[i]) + Sign(s[i - 1])) * std::min(1.0 * fabs(s[i - 1]), 0.5 * fabs(p[i]));
		}
		// Points in the middle
		else
		{
			p[i]  = (s[i - 1] * h[i] + s[i] * h[i - 1]) / (h[i - 1] + h[i]);
			dy[i] = (Sign(s[i - 1]) + Sign(s[i])) * std::min(1.0 * fabs(p[i]) / 2.0, std::min(1.0 * fabs(s[i]), 1.0 * fabs(s[i - 1])));
		}
	}

	// 3. a,b,c, and d
	for(unsigned int i = 0; i < N - 1; i++)
	{
		a.push_back((dy[i] + dy[i + 1] - 2.0 * s[i]) / pow(h[i], 2.0));
		b.push_back((3.0 * s[i] - 2.0 * dy[i] - dy[i + 1]) / h[i]);
		c.push_back(dy[i]);
		d.push_back(function_values[i]);
		if(std::isnan(a.back()) || std::isnan(b.back()) || std::isnan(c.back()) || std::isnan(d.back()))
			std::cout << "Warning in libphysica::Interpolation::Compute_Steffen_Coefficients(): Steffen coefficients in interpolation are NAN." << std::endl;
	}
}

unsigned int Interpolation::Bisection(double x, int jLeft, int jRight)
{
	while((jRight - jLeft) > 1)
	{
		int jm = (jRight + jLeft) >> 1;
		if(x >= x_values[jm])
			jLeft = jm;
		else
			jRight = jm;
	}
	return jLeft;
}

unsigned int Interpolation::Hunt(double x)
{

	// 1. Hunting phase returns jd,ju which bracket j
	int dj = 1;
	int jd;
	unsigned int ju;
	// Hunt up
	if(x > x_values[jLast])
	{
		jd = jLast;
		ju = jd + dj;
		while(x > x_values[ju])
		{
			jd = ju;
			ju += dj;
			// Check if we ran off the range:
			if(ju > N - 1)
			{
				ju = N - 1;
				break;
			}
			else
				dj += dj;
		}
	}
	// Hunt down
	else if(x < x_values[jLast])
	{
		ju = jLast;
		jd = ju - dj;
		while(x < x_values[jd])
		{
			ju = jd;
			jd -= dj;
			// Check if we ran off the range:
			if(jd < 0)
			{
				jd = 0;
				break;
			}
			else
				dj += dj;
		}
	}
	else
	{
		return jLast;
	}

	// 2. Bisection phase
	if((ju - jd) > 1)
		jd = Bisection(x, jd, ju);

	return jd;
}

// Find j such that list[j]<x<list[j+1]
unsigned int Interpolation::Locate(double x)
{
	unsigned int j;
	// Check if we are inside the domain, but allow a bit of extrapolation outside it.
	if(x < domain[0] || x > domain[1])
	{
		double boundary_tolerance_left	= 1e-2 * (x_values[1] - x_values[0]);
		double boundary_tolerance_right = 1e-2 * (x_values[N - 1] - x_values[N - 2]);
		if(fabs(x - domain[0]) < boundary_tolerance_left)
			j = 0;
		else if(fabs(x - domain[1]) < boundary_tolerance_right)
			j = N - 2;
		else
		{
			printf("\nError in libphysica::Interpolation::Locate(): x = %e lies outside the domain [%e,%e].\n\n", x, domain[0], domain[1]);
			std::exit(EXIT_FAILURE);
		}
	}
	else
	{
		// Use Bisection() or the Hunt method, depending of the last calls were correlated.
		j = correlated_calls ? Hunt(x) : Bisection(x, 0, N - 1);
	}
	// Check if the points are still correlated.
	correlated_calls = (fabs(j - jLast) < 10);
	jLast			 = j;
	return j;
}

// Constructors
Interpolation::Interpolation()
{
	std::vector<double> x_val = {-1.0, 0.0, 1.0};
	std::vector<double> y_val = {0.0, 0.0, 0.0};

	*this = Interpolation(x_val, y_val);
}

Interpolation::Interpolation(const std::vector<double>& arg_values, const std::vector<double>& func_values, double x_dim, double f_dim)
: N(arg_values.size()), x_values(arg_values), function_values(func_values), prefactor(1.0), jLast(0), correlated_calls(false)
{
	// Some initial checks
	if(x_values.size() != function_values.size())
	{
		std::cerr << "Error in libphysica::Interpolation::Interpolation(): Unequal length of argument and function lists: " << x_values.size() << " vs " << function_values.size() << std::endl;
		std::exit(EXIT_FAILURE);
	}
	for(unsigned int i = 1; i < N; i++)
	{
		if(x_values[i] <= x_values[i - 1])
		{
			std::cerr << "Error in libphysica::Interpolation::Interpolation(): Argument list not strictly increasing." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	// Transform units
	if(x_dim > 0.0)
		for(unsigned int i = 0; i < N; i++)
			x_values[i] *= x_dim;
	if(f_dim > 0.0)
		for(unsigned int i = 0; i < N; i++)
			function_values[i] *= f_dim;

	domain = {x_values[0], x_values[N - 1]};
	Compute_Steffen_Coefficients();
}

Interpolation::Interpolation(const std::vector<std::vector<double>>& data, double x_dim, double f_dim)
{
	std::vector<double> x;
	std::vector<double> f;
	for(unsigned int i = 0; i < data.size(); i++)
	{
		if(data[i].size() != 2)
		{
			std::cerr << "Error in libphysica::Interpolation::Interpolation(): Data table is has faulty dimensions (" << data[i].size() << ")" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		x.push_back(data[i][0]);
		f.push_back(data[i][1]);
	}
	*this = Interpolation(x, f, x_dim, f_dim);
}

void Interpolation::Set_Prefactor(double factor)
{
	prefactor = factor;
}

double Interpolation::Interpolate(double x)
{
	int j		= Locate(x);
	double x_j	= x_values[j];
	double inte = prefactor * (a[j] * pow((x - x_j), 3.0) + b[j] * pow((x - x_j), 2.0) + c[j] * (x - x_j) + d[j]);
	return inte;
}

void Interpolation::Multiply(double factor)
{
	prefactor *= factor;
}

// Function properties
double Interpolation::Derivative(double x, unsigned int derivation)
{
	int j	   = Locate(x);
	double x_j = x_values[j];
	if(derivation == 0)
		return Interpolate(x);
	else if(derivation == 1)
		return prefactor * (3.0 * a[j] * pow((x - x_j), 2.0) + 2.0 * b[j] * (x - x_j) + c[j]);
	else if(derivation == 2)
		return prefactor * (6.0 * a[j] * (x - x_j) + 2.0 * b[j]);
	else if(derivation == 3)
		return prefactor * (6.0 * a[j]);
	else
		return 0.0;
}

double Interpolation::Integrate(double x_1, double x_2)
{
	double sign = +1.0;
	if(x_1 > x_2)
	{
		std::swap(x_1, x_2);
		sign = -1.0;
	}
	int i_1			= Locate(x_1);
	int i_2			= Locate(x_2);
	double integral = 0;
	for(int i = 0; i < (i_2 - i_1 + 1); i++)
	{
		int j				  = i_1 + i;
		double x_j			  = x_values[j];
		double x_left		  = (i == 0) ? x_1 : x_j;
		double x_right		  = (i == (i_2 - i_1)) ? x_2 : x_values[j + 1];
		double stemfunc_left  = prefactor * (a[j] / 4.0 * pow((x_left - x_j), 4.0) + b[j] / 3.0 * pow((x_left - x_j), 3.0) + c[j] / 2.0 * pow((x_left - x_j), 2.0) + d[j] * x_left);
		double stemfunc_right = prefactor * (a[j] / 4.0 * pow((x_right - x_j), 4.0) + b[j] / 3.0 * pow((x_right - x_j), 3.0) + c[j] / 2.0 * pow((x_right - x_j), 2.0) + d[j] * x_right);
		integral += stemfunc_right - stemfunc_left;
	}
	return sign * integral;
}

double Interpolation::Local_Minimum(double x_1, double x_2)
{
	libphysica::Check_For_Error(x_2 < x_1, "Interpolation::Local_Minimum()", "Faulty order of arguments.");
	double f_left  = Interpolate(x_1);
	double f_right = Interpolate(x_2);
	int i_1		   = Locate(x_1);
	int i_2		   = Locate(x_2);
	if(i_1 == i_2)
		return std::min(f_left, f_right);
	else
	{
		// Find the smallest value of function_values between i_1+1 and i_2.
		double min_entry = *std::min_element(function_values.begin() + i_1 + 1, function_values.begin() + i_2);
		return std::min({f_left, min_entry, f_right});
	}
}

double Interpolation::Local_Maximum(double x_1, double x_2)
{
	libphysica::Check_For_Error(x_2 < x_1, "Interpolation::Local_Minimum()", "Faulty order of arguments.");
	double f_left  = Interpolate(x_1);
	double f_right = Interpolate(x_2);
	int i_1		   = Locate(x_1);
	int i_2		   = Locate(x_2);
	if(i_1 == i_2)
		return std::max(f_left, f_right);
	else
	{
		// Find the largest value of function_values between i_1+1 and i_2.
		double max_entry = *std::max_element(function_values.begin() + i_1 + 1, function_values.begin() + i_2);
		return std::max({f_left, max_entry, f_right});
	}
}

double Interpolation::Global_Minimum()
{
	return *std::min_element(function_values.begin(), function_values.end());
}

double Interpolation::Global_Maximum()
{
	return *std::max_element(function_values.begin(), function_values.end());
}

void Interpolation::Save_Function(std::string filename, unsigned int points)
{
	std::ofstream f;
	f.open(filename);
	std::vector<double> x_points = Linear_Space(domain[0], domain[1], points);
	for(auto& x : x_points)
		f << x << "\t" << Interpolate(x) << std::endl;
	f.close();
}

// 1.2 Two-dimensional interpolation (bilinear interpolation)
Interpolation_2D::Interpolation_2D()
{
	std::vector<double> x_val = {-1.0, 0.0, 1.0};
	std::vector<double> y_val = {-1.0, 0.0, 1.0};
	std::vector<std::vector<double>> f_val(3, std::vector<double>(3, 0.0));

	*this = Interpolation_2D(x_val, y_val, f_val);
}

Interpolation_2D::Interpolation_2D(std::vector<double> x_val, std::vector<double> y_val, std::vector<std::vector<double>> func_values, double x_dim, double y_dim, double f_dim)
: N_x(x_val.size()), N_y(y_val.size()), x_values(x_val), y_values(y_val), function_values(func_values), prefactor(1.0)
{
	// Transform units
	if(x_dim > 0.0)
		for(unsigned int i = 0; i < N_x; i++)
			x_values[i] *= x_dim;
	if(y_dim > 0.0)
		for(unsigned int i = 0; i < N_y; i++)
			y_values[i] *= y_dim;
	if(f_dim > 0.0)
		for(unsigned int i = 0; i < N_x; i++)
			for(unsigned int j = 0; j < N_y; j++)
				function_values[i][j] *= f_dim;

	// Create dummy 1D interpolations to use for index location.
	std::vector<double> dummy_x(N_x, 0.0);
	std::vector<double> dummy_y(N_y, 0.0);
	x_int = Interpolation(x_values, std::vector<double>(N_x, 0.0));
	y_int = Interpolation(y_values, std::vector<double>(N_y, 0.0));

	domain = {x_int.domain, y_int.domain};
}

Interpolation_2D::Interpolation_2D(std::vector<std::vector<double>> data_table, double x_dim, double y_dim, double f_dim)
{
	std::vector<double> x;
	std::vector<double> y;

	for(unsigned int i = 0; i < data_table.size(); i++)
	{
		if(data_table[i].size() != 3)
		{
			std::cerr << "Error in libphysica::Interpolation_2D::Interpolation_2D(): Data table is has faulty dimensions (" << data_table[i].size() << ")" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		x.push_back(data_table[i][0]);
		y.push_back(data_table[i][1]);
	}

	// Delete dublicates in x and y list and sort
	std::sort(x.begin(), x.end());
	std::sort(y.begin(), y.end());
	x.erase(unique(x.begin(), x.end()), x.end());
	y.erase(unique(y.begin(), y.end()), y.end());

	if(x.size() * y.size() != data_table.size())
	{
		std::cerr << "Error in libphysica::Interpolation_2D::Interpolation_2D(): List lenghts do not fit." << std::endl;
		std::exit(EXIT_FAILURE);
	}

	std::vector<std::vector<double>> f(x.size(), std::vector<double>(y.size()));
	unsigned int i = 0;
	for(unsigned int i_x = 0; i_x < x.size(); i_x++)
		for(unsigned int i_y = 0; i_y < y.size(); i_y++)
		{
			if(x[i_x] != data_table[i][0] || y[i_y] != data_table[i][1])
			{
				std::cerr << "Error in libphysica::Interpolation_2D::Interpolation_2D(): Data table was not in right format." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			f[i_x][i_y] = data_table[i][2];
			i++;
		}

	*this = Interpolation_2D(x, y, f, x_dim, y_dim, f_dim);
}

double Interpolation_2D::Interpolate(double x, double y)
{
	unsigned int i = x_int.Locate(x);
	unsigned int j = y_int.Locate(y);

	double t = (x - x_values[i]) / (x_values[i + 1] - x_values[i]);
	double u = (y - y_values[j]) / (y_values[j + 1] - y_values[j]);

	double f0 = function_values[i][j];
	double f1 = function_values[i + 1][j];
	double f2 = function_values[i + 1][j + 1];
	double f3 = function_values[i][j + 1];

	return (1.0 - t) * (1.0 - u) * f0 + t * (1.0 - u) * f1 + t * u * f2 + (1.0 - t) * u * f3;
}

void Interpolation_2D::Set_Prefactor(double factor)
{
	prefactor = factor;
}

void Interpolation_2D::Multiply(double factor)
{
	prefactor *= factor;
}

// Function properties
double Interpolation_2D::Global_Minimum()
{
	std::vector<double> row_minima;
	for(auto& row : function_values)
		row_minima.push_back(*std::min_element(row.begin(), row.end()));
	return *std::min_element(row_minima.begin(), row_minima.end());
}
double Interpolation_2D::Global_Maximum()
{
	std::vector<double> row_maxima;
	for(auto& row : function_values)
		row_maxima.push_back(*std::max_element(row.begin(), row.end()));
	return *std::max_element(row_maxima.begin(), row_maxima.end());
}

void Interpolation_2D::Save_Function(std::string filename, unsigned int x_points, unsigned int y_points)
{
	std::ofstream f;
	f.open(filename);
	if(y_points == 0)
		y_points = x_points;
	std::vector<double> x_list = Linear_Space(domain[0][0], domain[0][1], x_points);
	std::vector<double> y_list = Linear_Space(domain[1][0], domain[1][1], y_points);
	for(auto& x : x_list)
		for(auto& y : y_list)
			f << x << "\t" << y << "\t" << Interpolate(x, y) << std::endl;
	f.close();
}

// 2. Root finding
// Root finding with Ridder's method
double Find_Root(std::function<double(double)> func, double xLeft, double xRight, double xAccuracy)
{
	const int Max_Iterations = 50;
	// 1. Check if xLeft<xRight, otherwise swap.
	if(xLeft > xRight)
	{
		double temp = xLeft;
		xLeft		= xRight;
		xRight		= temp;
	}

	// 2. Compute functions at boundary
	double fLeft  = func(xLeft);
	double fRight = func(xRight);

	// 3. Check if xLeft and xRight bracket a root or already yield a root. Also check for NaN's.
	if(std::isnan(fLeft) || std::isnan(fRight))
	{
		std::cerr << "Error in libphysica::Find_Root(): Function returns nan at the brackets." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(fLeft * fRight >= 0.0)
	{
		if(fLeft == 0)
			return xLeft;
		else if(fRight == 0)
			return xRight;
		else
		{
			std::cerr << "Error in libphysica::Find_Root(): f(xLeft = " << xLeft << ") * f(xRight = " << xRight << ") = (" << fLeft << ") * (" << fRight << ") > 0.0" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	// 4. Ridder's method
	else
	{
		double x1	  = xLeft;
		double x2	  = xRight;
		double f1	  = fLeft;
		double f2	  = fRight;
		double result = -9.9e99;
		for(int i = 0; i < Max_Iterations; i++)
		{
			// Mid point
			double x3 = (x1 + x2) / 2.0;

			double f3 = func(x3);
			// New point
			double x4 = x3 + (x3 - x1) * Sign(f1 - f2) * f3 / sqrt(f3 * f3 - f1 * f2);
			// Check if we found the root
			if(fabs(x4 - result) < xAccuracy)
				return x4;
			// Prepare next iteration
			result	  = x4;
			double f4 = func(x4);
			if(f4 == 0.0)
				return result;
			// a) x3 and x4 bracket the root
			if(Sign(f3, f4) != f3)
			{
				x1 = x3;
				f1 = f3;
				x2 = x4;
				f2 = f4;
			}
			// b) x1 and x4 bracket the root
			else if(Sign(f1, f4) != f1)
			{
				x2 = x4;
				f2 = f4;
			}
			// c) x2 and x4 bracket the root
			else if(Sign(f2, f4) != f2)
			{
				x1 = x4;
				f1 = f4;
			}
			else
			{
				std::cerr << "Error in libphysica::Find_Root(). Ridder's method does not reach the root." << std::endl;
				std::exit(EXIT_FAILURE);
			}
		}
		std::cout << "Warning in libphysica::Find_Root(): Iterations exceed the maximum. Final value f(" << result << ")=" << func(result) << std::endl;
		return result;
	}
}

// 3. Minimization
// 3.1 One-dimensional
// Find a<b<c such that f(b)<f(a) and f(b) < f(c).
struct Bracket_Method
{
	double ax, bx, cx, fa, fb, fc;

	template <class T>
	void Bracket(const double a, const double b, T& func)
	{
		const double golden_ratio = 1.618034;
		double GLIMIT			  = 100.0;
		double TINY				  = 1.0e-20;

		ax = a;
		bx = b;
		double fu;
		fa = func(ax);
		fb = func(bx);
		if(fb > fa)
		{
			std::swap(ax, bx);
			std::swap(fb, fa);
		}
		cx = bx + golden_ratio * (bx - ax);
		fc = func(cx);
		while(fb > fc)
		{
			double r	= (bx - ax) * (fb - fc);
			double q	= (bx - cx) * (fb - fa);
			double u	= bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * Sign(std::max(std::fabs(q - r), TINY), q - r));
			double ulim = bx + GLIMIT * (cx - bx);
			if((bx - u) * (u - cx) > 0.0)
			{
				fu = func(u);
				if(fu < fc)
				{
					ax = bx;
					bx = u;
					fa = fb;
					fb = fu;
					return;
				}
				else if(fu > fb)
				{
					cx = u;
					fc = fu;
					return;
				}
				u  = cx + golden_ratio * (cx - bx);
				fu = func(u);
			}
			else if((cx - u) * (u - ulim) > 0.0)
			{
				fu = func(u);
				if(fu < fc)
				{
					Shift3(bx, cx, u, u + golden_ratio * (u - bx));
					Shift3(fb, fc, fu, func(u));
				}
			}
			else if((u - ulim) * (ulim - cx) >= 0.0)
			{
				u  = ulim;
				fu = func(u);
			}
			else
			{
				u  = cx + golden_ratio * (cx - bx);
				fu = func(u);
			}
			Shift3(ax, bx, cx, u);
			Shift3(fa, fb, fc, fu);
		}
	}

	inline void Shift2(double& a, double& b, const double c)
	{
		a = b;
		b = c;
	}
	inline void Shift3(double& a, double& b, double& c, const double d)
	{
		a = b;
		b = c;
		c = d;
	}

	inline void Move3(double& a, double& b, double& c, const double d, const double e, const double f)
	{
		a = d;
		b = e;
		c = f;
	}
};

struct Brent : Bracket_Method
{
	double x_min, f_min;
	const double tol;

	Brent(const double toll = 3e-8)
	: tol(toll) {};

	template <class T>
	double Minimize(T& func)
	{
		const int ITMAX	   = 100;
		const double CGOLD = 0.3819660;
		const double ZEPS  = std::numeric_limits<double>::epsilon();
		double a, b, d = 0.0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm, e = 0.0;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = func(x);
		for(int iter = 0; iter < ITMAX; iter++)
		{
			xm	 = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = tol * std::fabs(x) + ZEPS);
			if(std::fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
			{
				f_min = fx;
				x_min = x;
				return x_min;
			}
			if(std::fabs(e) > tol1)
			{
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if(q > 0.0)
					p = -p;
				q	  = std::fabs(q);
				etemp = e;
				e	  = d;
				if(std::fabs(p) >= std::fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				else
				{
					d = p / q;
					u = x + d;
					if(u - a < tol2 || b - u < tol2)
						d = Sign(tol1, xm - x);
				}
			}
			else
				d = CGOLD * (e = (x >= xm ? a - x : b - x));

			u  = (std::fabs(d) >= tol1 ? x + d : x + Sign(tol1, d));
			fu = func(u);

			if(fu <= fx)
			{
				if(u >= x)
					a = x;
				else
					b = x;
				Shift3(v, w, x, u);
				Shift3(fv, fw, fx, fu);
			}
			else
			{
				if(u < x)
					a = u;
				else
					b = u;
				if(fu <= fw || w == x)
				{
					v  = w;
					w  = u;
					fv = fw;
					fw = fu;
				}
				else if(fu <= fv || v == x || v == w)
				{
					v  = u;
					fv = fu;
				}
			}
		}
		std::cerr << "Error in Brent::Minimize(): Too many iterations." << std::endl;
		std::exit(EXIT_FAILURE);
	}
};

double Find_Maximum(std::function<double(double)> func, double xLeft, double xRight, double tol)
{
	auto minus_func = [&func](double x) {
		return -1.0 * func(x);
	};
	return Find_Minimum(minus_func, xLeft, xRight, tol);
}

double Find_Minimum(std::function<double(double)> func, double xLeft, double xRight, double tol)
{
	Brent brent(tol);
	brent.Bracket(xLeft, xRight, func);
	double x_min = brent.Minimize(func);
	return x_min;
}

// 3.2 Multi-dimensional
std::vector<double> Minimization::minimize(std::vector<double>& starting_point, const double delta, std::function<double(std::vector<double>)> func)
{
	std::vector<double> deltas(starting_point.size(), delta);
	return minimize(starting_point, deltas, func);
}

std::vector<double> Minimization::minimize(std::vector<double>& starting_point, std::vector<double>& deltas, std::function<double(std::vector<double>)> func)
{
	int ndim = starting_point.size();
	std::vector<std::vector<double>> pp(ndim + 1, std::vector<double>(ndim, 0.0));
	for(int i = 0; i < ndim + 1; i++)
	{
		for(int j = 0; j < ndim; j++)
			pp[i][j] = starting_point[j];
		if(i != 0)
			pp[i][i - 1] += deltas[i - 1];
	}
	return minimize(pp, func);
}

std::vector<double> Minimization::minimize(std::vector<std::vector<double>>& pp, std::function<double(std::vector<double>)> func)
{
	const int NMAX	  = 5000;
	const double TINY = 1.0e-10;
	mpts			  = pp.size();		// rows
	ndim			  = pp[0].size();	// columns
	std::vector<double> psum(ndim), pmin(ndim), x(ndim);
	current_simplex = pp;
	y.resize(mpts);
	for(int i = 0; i < mpts; i++)
	{
		for(int j = 0; j < ndim; j++)
			x[j] = current_simplex[i][j];
		y[i] = func(x);
	}

	nfunc = 0;
	get_psum(current_simplex, psum);
	for(;;)
	{
		int ihi, ilo, inhi;
		ilo = 0;
		// First we must determine which point is the highest (worst), next-highest, and lowest (best), by looping over the points in the simplex.
		ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
		for(int i = 0; i < mpts; i++)
		{
			if(y[i] <= y[ilo])
				ilo = i;
			if(y[i] > y[ihi])
			{
				inhi = ihi;
				ihi	 = i;
			}
			else if(y[i] > y[inhi] && i != ihi)
				inhi = i;
		}
		double rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]) + TINY);
		// Compute the fractional range from highest to lowest and return if satisfactory.
		if(rtol < ftol)
		{
			std::swap(y[0], y[ilo]);
			for(int i = 0; i < ndim; i++)
			{
				std::swap(current_simplex[0][i], current_simplex[ilo][i]);
				pmin[i] = current_simplex[0][i];
			}
			fmin = y[0];
			return pmin;
		}
		if(nfunc >= NMAX)
		{
			std::cerr << "Error in libphysica::Minimization::minimize(): NMAX exceeded." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		nfunc += 2;

		// Begin a new iteration. First extrapolate by a factor 􏰱1 through the face of the simplex across from the high point, i.e., reflect the simplex from the high point.
		double ytry = amotry(current_simplex, y, psum, ihi, -1.0, func);
		if(ytry <= y[ilo])
			ytry = amotry(current_simplex, y, psum, ihi, 2.0, func);   // Gives a result better than the best point, so try an additional extrapolation by a factor 2.
		else if(ytry >= y[inhi])
		{
			// The reflected point is worse than the second-highest, so look for an interme- diate lower point, i.e., do a one-dimensional contraction.
			double ysave = y[ihi];
			ytry		 = amotry(current_simplex, y, psum, ihi, 0.5, func);
			if(ytry >= ysave)
			{
				// Can’t seem to get rid of that high point.
				for(int i = 0; i < mpts; i++)
				{
					if(i != ilo)
					{
						for(int j = 0; j < ndim; j++)
							current_simplex[i][j] = psum[j] = 0.5 * (current_simplex[i][j] + current_simplex[ilo][j]);
						y[i] = func(psum);
					}
				}
				nfunc += ndim;					   // Keep track of function evaluations.
				get_psum(current_simplex, psum);   // Recompute psum.
			}
		}
		else
			--nfunc;   // Correct the evaluation count.
	}				   // Go back to the test of doneness and the next iteration
}

void Minimization::get_psum(std::vector<std::vector<double>>& p, std::vector<double>& psum)	  // Utility function.
{
	for(int j = 0; j < ndim; j++)
	{
		double sum = 0.0;
		for(int i = 0; i < mpts; i++)
			sum += p[i][j];
		psum[j] = sum;
	}
}

double Minimization::amotry(std::vector<std::vector<double>>& p, std::vector<double>& y, std::vector<double>& psum, const int ihi, const double fac, std::function<double(std::vector<double>)> func)
{
	std::vector<double> ptry(ndim);
	double fac1 = (1.0 - fac) / ndim;
	double fac2 = fac1 - fac;
	for(int j = 0; j < ndim; j++)
		ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
	double ytry = func(ptry);	// Evaluate the function at the trial point.
	if(ytry < y[ihi])			// If it’s better than the highest, then replace the highest.
	{
		y[ihi] = ytry;
		for(int j = 0; j < ndim; j++)
		{
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	return ytry;
}

}	// namespace libphysica
