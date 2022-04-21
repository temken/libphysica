// Disclaimer:
// Some of the function implementations were made with the help of the
//"Numerical Recipes 3rd Edition: The Art of Scientific Computing"
// by William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery

#include "libphysica/Special_Functions.hpp"

#include <cmath>
#include <limits>	// std::numeric_limits

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "libphysica/Integration.hpp"
#include "libphysica/Numerics.hpp"

namespace libphysica
{
using namespace std::complex_literals;

// 1. Simple functions
int Sign(double arg)
{
	if(arg > 0.0)
		return 1;
	else if(arg == 0.0)
		return 0;
	else
		return -1;
}

double Sign(double x, double y)
{
	if(Sign(x) == Sign(y))
		return x;
	else
		return -1.0 * x;
}

double StepFunction(double x)
{
	if(x >= 0)
		return 1.0;
	else
		return 0.0;
}

double Round(double N, unsigned int digits)
{
	if(N == 0)
		return 0;
	unsigned int digits_max = 7;
	if(digits > digits_max)
	{
		std::cerr << "Error in libphysica::Round(): Significant digits > " << digits_max << "." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Make the argument a positive number.
	double sign = Sign(N);
	N *= sign;

	// Cut off the decimal power
	double DecimalPower = floor(log10(N));
	double prefactor	= N * pow(10, -DecimalPower);

	// Round the prefactor
	prefactor = std::floor(prefactor * pow(10.0, digits - 1) + 0.5);
	prefactor = prefactor * pow(10.0, -1.0 * digits + 1);

	return sign * prefactor * pow(10, DecimalPower);
}

Vector Round(const Vector& vec, unsigned int digits)
{
	Vector rounded_vec = vec;
	for(unsigned int i = 0; i < vec.Size(); i++)
		rounded_vec[i] = Round(vec[i], digits);
	return rounded_vec;
}

Matrix Round(const Matrix& matrix, unsigned int digits)
{
	Matrix rounded_matrix = matrix;
	for(unsigned int i = 0; i < matrix.Rows(); i++)
		for(unsigned int j = 0; j < matrix.Columns(); j++)
			rounded_matrix[i][j] = Round(matrix[i][j], digits);
	return rounded_matrix;
}

double Relative_Difference(double a, double b)
{
	double d   = std::fabs(a - b);
	double max = std::max(fabs(a), fabs(b));
	return d / max;
}

bool Floats_Equal(double a, double b, double tol)
{
	if(Relative_Difference(a, b) < tol)
		return true;
	else
		return false;
}

// 2. Special functions
// 2.1 Gamma functions
std::vector<double> FactorialList = {1.0};
double Factorial(unsigned int n)
{
	if(n > 170)
	{
		std::cerr << "Error in libphysica::Factorial: Overflow for " << n << "!." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(n < FactorialList.size())
		return FactorialList[n];
	else
	{
		while(FactorialList.size() <= n)
			FactorialList.push_back(FactorialList.back() * FactorialList.size());
		return FactorialList.back();
	}
}

double Binomial_Coefficient(int n, int k)
{
	if(k < 0 || n < 0)
	{
		std::cerr << "Warning in libphysica::BinomialCoefficient(): negative arguments. Return 0." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(n < k)
		return 0;
	else if(n > 170)
	{
		return floor(0.5 + exp(GammaLn(n + 1.0) - GammaLn(k + 1.0) - GammaLn(n - k + 1.0)));
	}
	else
		return floor(0.5 + Factorial(n) / Factorial(k) / Factorial(n - k));
}

// Logarithmic gamma function
double GammaLn(double x)
{
	double cof[14] = {57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3, .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};
	if(x <= 0)
	{
		std::cerr << "Error in libphysica::GammaLn(x): x<=0." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	double sum = 0.999999999999997092;
	double y   = x;
	double tmp = x + 671.0 / 128.0;
	tmp		   = (x + 0.5) * log(tmp) - tmp;
	for(int j = 0; j < 14; j++)
	{
		sum += cof[j] / ++y;
	}
	return tmp + log(2.5066282746310005 * sum / x);
}

double Gamma(double x)
{
	return exp(GammaLn(x));
}

double Upper_Incomplete_Gamma(double x, double s)
{
	return Gamma(s) * GammaQ(x, s);
}

double Lower_Incomplete_Gamma(double x, double s)
{
	return Gamma(s) * GammaP(x, s);
}

// Q(x,a) via integration;
double GammaQint(double x, double a)
{
	// Compute P(x,a)
	double gammaP;
	double gln = GammaLn(a);
	// How far to integrate N sqrt(a) around the peak at a-1:
	double N	 = 10;
	double tPeak = a - 1.0;
	double tMin	 = std::max(0.0, tPeak - N * sqrt(a));
	double tMax	 = tPeak + N * sqrt(a);
	// If x lies far away from the peak:
	if(x > tMax)
		gammaP = 1.0;
	else if(x < tMin)
		gammaP = 0.0;
	// Numerical integration
	else
	{

		// integrand
		std::function<double(double)> integrand = [a, gln](double t) {
			return exp(-gln - t + log(t) * (a - 1.0));
		};
		if(x < tMin)
			tMin = 0.0;
		// Precision
		double eps = Find_Epsilon(integrand, tMin, x, 1e-5);
		// Integrate
		gammaP = Integrate(integrand, tMin, x, eps);
	}

	return 1.0 - gammaP;
}

// Series expansion of P(x,a)
double GammaPser(double x, double a)
{
	double eps = std::numeric_limits<double>::epsilon();
	double sum, del, ap;
	double gln = GammaLn(a);
	ap		   = a;
	del = sum = 1.0 / a;
	while(fabs(del) > fabs(sum) * eps)
	{
		ap++;
		del *= x / ap;
		sum += del;
	}
	return sum * exp(-x + a * log(x) - gln);
}

// Continued fraction representation of Q(x,a)
double GammaQcf(double x, double a)
{
	// Precision
	double eps	 = std::numeric_limits<double>::epsilon();
	double FPMIN = std::numeric_limits<double>::min() / eps;
	double del	 = 0.0;
	double gln	 = GammaLn(a);
	double b	 = x + 1.0 - a;
	double c	 = 1.0 / FPMIN;
	double d	 = 1.0 / b;
	double h	 = d;
	int i		 = 1;
	while(fabs(del - 1.0) > eps)
	{
		double an = -1.0 * i * (i - a);
		b += 2.0;
		d = an * d + b;
		if(fabs(d) < FPMIN)
			d = FPMIN;
		c = b + an / c;
		if(fabs(c) < FPMIN)
			c = FPMIN;
		d	= 1.0 / d;
		del = d * c;
		h *= del;
	}
	return exp(-x + a * log(x) - gln) * h;
}

// Final function using different methods for different parts of the domain
double GammaQ(double x, double a)
{
	double aMax = 100.0;
	if(x < 0.0 || a <= 0.0)
	{
		std::cerr << "Error in libphysica::GammaQ(" << x << "," << a << "): Invalid arguments." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(x == 0)
		return 1.0;
	else if(a > aMax)
		return GammaQint(x, a);
	else if(x < a + 1.0)
		return 1.0 - GammaPser(x, a);
	else
		return GammaQcf(x, a);
}

double GammaP(double x, double a)
{
	return 1.0 - GammaQ(x, a);
}

// Inverse incomplete gamma function. (Solves P(x,a)=p for x.)
double Inv_GammaP(double p, double a)
{
	// Check the arguments
	if(a <= 0.0)
	{
		std::cerr << "Error in libphysica::Inv_GammaP(): a must be positive." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if(p >= 1.0)
		return std::max(100.0, a + 100. * sqrt(a));
	if(p <= 0.0)
		return 0.0;

	// Parameter
	double x;
	double gln	= GammaLn(a);
	double a1	= a - 1.0;
	double lna1 = log(a1);
	double afac = exp(a1 * (lna1 - 1.0) - gln);

	// Initial guess 1
	if(a > 1.0)
	{
		double pp = (p < 0.5) ? p : 1.0 - p;
		double t  = sqrt(-2.0 * log(pp));
		x		  = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
		if(p < 0.5)
			x = -x;
		x = std::max(1.0e-3, a * pow(1.0 - 1.0 / (9. * a) - x / (3. * sqrt(a)), 3.0));
	}
	// Initial guess 2
	else
	{
		double t = 1.0 - a * (0.253 + a * 0.12);
		if(p < t)
			x = pow(p / t, 1. / a);
		else
			x = 1. - log(1. - (p - t) / (1. - t));
	}
	// Halley's method
	double EPS = 1.0e-8;
	for(int i = 0; i < 12; i++)
	{
		if(x <= 0.0)
			return 0.0;
		double error = GammaP(x, a) - p;
		double t;
		if(a > 1.0)
			t = afac * exp(-(x - a1) + a1 * (log(x) - lna1));
		else
			t = exp(-x + a1 * log(x) - gln);
		double u = error / t;
		x -= (t = u / (1. - 0.5 * std::min(1., u * ((a - 1.) / x - 1))));
		if(x <= 0.)
			x = 0.5 * (x + t);
		if(fabs(t) < EPS * x)
			break;
	}
	return x;
}

// Solves Q(x,a)=p for x.
double Inv_GammaQ(double q, double a)
{
	return Inv_GammaP(1.0 - q, a);
}

// 2.2 Scalar spherical Harmonics
std::complex<double> Spherical_Harmonics(int l, int m, double theta, double phi)
{
	return boost::math::spherical_harmonic(l, m, theta, phi);
}

// 2.3 Vector spherical Harmonics
std::complex<double> VSH_Y_Component(int component, int l, int m, int l_hat, int m_hat)
{
	switch(component)
	{
		case 0:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 && m_hat == m + 1)
				return -1.0 / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 && m_hat == m - 1)
				return 1.0 / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m + 1)
				return 1.0 / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m - 1)
				return -1.0 / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 1:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 && m_hat == m + 1)
				return 1.0i / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 && m_hat == m - 1)
				return 1.0i / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 && m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 2:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m))
				return 0.0;
			else if(l_hat == l + 1)
				return 1.0 * std::sqrt(1.0 * (l - m + 1) * (l + m + 1) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1)
				return 1.0 * std::sqrt(1.0 * (l - m) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		default:
			std::cerr << "Error in VSH_Y_Component(): Component " << component << " out of bound." << std::endl;
			std::exit(EXIT_FAILURE);
	}
	return 0.0;
}

std::vector<std::complex<double>> Vector_Spherical_Harmonics_Y(int l, int m, double theta, double phi)
{
	std::vector<std::complex<double>> Y(3, 0.0);
	for(unsigned int i = 0; i < 3; i++)
		for(int l_hat = l - 1; l_hat < l + 2; l_hat += 2)
			for(int m_hat = m - 1; m_hat < m + 2; m_hat++)
				if(std::abs(m_hat) <= l_hat)
					Y[i] += VSH_Y_Component(i, l, m, l_hat, m_hat) * Spherical_Harmonics(l_hat, m_hat, theta, phi);
	return Y;
}

std::complex<double> VSH_Psi_Component(int component, int l, int m, int l_hat, int m_hat)
{
	switch(component)
	{
		case 0:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 and m_hat == m + 1)
				return l / 2.0 * std::sqrt(1.0 * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 and m_hat == m - 1)
				return -l / 2.0 * std::sqrt(1.0 * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m + 1)
				return (l + 1) / 2.0 * std::sqrt(1.0 * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m - 1)
				return -(l + 1) / 2.0 * std::sqrt(1.0 * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 1:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m - 1 && m_hat != m + 1))
				return 0.0;
			else if(l_hat == l + 1 and m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * l * l * (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l + 1 and m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * l * l * (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m + 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + 1) * (l + 1) * (l - m - 1) * (l - m) / (2 * l - 1) / (2 * l + 1));
			else if(l_hat == l - 1 and m_hat == m - 1)
				return -1.0i / 2.0 * std::sqrt(1.0 * (l + 1) * (l + 1) * (l + m - 1) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		case 2:
			if((l_hat != l - 1 && l_hat != l + 1) || (m_hat != m))
				return 0.0;
			else if(l_hat == l + 1)
				return -l * std::sqrt(1.0 * (l - m + 1) * (l + m + 1) / (2 * l + 3) / (2 * l + 1));
			else if(l_hat == l - 1)
				return (1 + l) * std::sqrt(1.0 * (l - m) * (l + m) / (2 * l - 1) / (2 * l + 1));
			break;
		default:
			std::cerr << "Error in VSH_Psi_Component(): Component " << component << " out of bound." << std::endl;
			std::exit(EXIT_FAILURE);
	}
	return 0.0;
}

std::vector<std::complex<double>> Vector_Spherical_Harmonics_Psi(int l, int m, double theta, double phi)
{
	std::vector<std::complex<double>> Psi(3, 0.0);
	for(unsigned int i = 0; i < 3; i++)
		for(int l_hat = l - 1; l_hat < l + 2; l_hat += 2)
			for(int m_hat = m - 1; m_hat < m + 2; m_hat++)
				if(std::abs(m_hat) <= l_hat)
					Psi[i] += VSH_Psi_Component(i, l, m, l_hat, m_hat) * Spherical_Harmonics(l_hat, m_hat, theta, phi);
	return Psi;
}

// 2.4 Other special functions
double Dawson_Integral(double x)
{
	double ans;
	static const double H = 0.4, A1 = 2.0 / 3.0, A2 = 0.4, A3 = 2.0 / 7.0;

	if(std::fabs(x) < 0.2)
	{
		double x2 = x * x;
		ans		  = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
	}
	else
	{
		static const int NMAX = 6;
		static std::vector<double> c(NMAX);
		for(int i = 0; i < NMAX; i++)
			c[i] = exp(-(2.0 * i + 1.0) * (2.0 * i + 1.0) * H * H);

		double xx  = std::fabs(x);
		int n0	   = 2 * int(0.5 * xx / H + 0.5);
		double xp  = xx - n0 * H;
		double e1  = std::exp(2.0 * xp * H);
		double e2  = e1 * e1;
		double d1  = n0 + 1;
		double d2  = d1 - 2.0;
		double sum = 0.0;
		for(int i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
			sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
		ans = 0.5641895835 * libphysica::Sign(std::exp(-xp * xp), x) * sum;
	}
	return ans;
}

double Erfi(double x)
{
	return 2.0 / std::sqrt(M_PI) * std::exp(x * x) * Dawson_Integral(x);
}

double Inv_Erf(double p)
{
	// return inverfc(1.-p);
	if(fabs(p - 1.0) < 1e-16)
	{
		std::cerr << "Warning in libphysica::Inv_erf(double): The argument p = " << p << " is very close to 1.0. Return 10." << std::endl;
		return 10.0;
	}
	else if(fabs(p) >= 1.0)
	{
		std::cerr << "Error in libphysica::Inv_erf(): Invalid argument |p| = |" << p << "| > 1" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::function<double(double)> fct = [p](double x) {
			return erf(x) - p;
		};
		double xLeft  = -10.0;
		double xRight = 10.0;
		return Find_Root(fct, xLeft, xRight, 1.0e-4);
	}
}
}	// namespace libphysica