#ifndef __Special_Functions_hpp_
#define __Special_Functions_hpp_

#include <complex>

#include "libphysica/Linear_Algebra.hpp"

namespace libphysica
{
// 1. Simple functions
extern int Sign(double arg);
extern double Sign(double x, double y);	  // Returns x with the sign of y.
extern double StepFunction(double x);
extern double Round(double N, unsigned int digits = 3);
extern Vector Round(const Vector& vec, unsigned int digits = 3);
extern Matrix Round(const Matrix& matrix, unsigned int digits = 3);
extern double Relative_Difference(double a, double b);
extern bool Floats_Equal(double a, double b, double tol = 1e-10);

// 2. Special functions
// 2.1 Gamma functions
extern double Factorial(unsigned int n);
extern double Binomial_Coefficient(int n, int k);
extern double GammaLn(double x);
extern double Gamma(double x);
extern double Upper_Incomplete_Gamma(double x, double s);
extern double Lower_Incomplete_Gamma(double x, double s);
extern double GammaQ(double x, double a);		// Incomplete Gamma Q(x,a)
extern double GammaP(double x, double a);		// Incomplete Gamma P(x,a)
extern double Inv_GammaP(double p, double a);	// Solves P(x,a)=p for x.
extern double Inv_GammaQ(double q, double a);	// Solves Q(x,a)=q for x.

// 2.2 Scalar spherical Harmonics
extern std::complex<double> Spherical_Harmonics(int l, int m, double theta, double phi);

// 2.3 Vector spherical Harmonics
extern std::complex<double> VSH_Y_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Y(int l, int m, double theta, double phi);

extern std::complex<double> VSH_Psi_Component(int component, int l, int m, int l_hat, int m_hat);
extern std::vector<std::complex<double>> Vector_Spherical_Harmonics_Psi(int l, int m, double theta, double phi);

// 2.4 Other special functions
extern double Erfi(double x);
extern double Inv_Erf(double p);

}	// namespace libphysica

#endif
