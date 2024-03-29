#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <random>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"

using namespace libphysica;
using namespace libphysica::natural_units;
using namespace std::complex_literals;

// 1. Simple functions
TEST(TestSpecialFunctions, TestSign)
{
	// ARRANGE
	double positive_number = 10.3;
	double negative_number = -2.1;
	// ACT & ASSERT
	ASSERT_EQ(Sign(positive_number), 1);
	ASSERT_EQ(Sign(negative_number), -1);
}

TEST(TestSpecialFunctions, TestSign2)
{
	// ARRANGE
	double positive_number = 10.3;
	double negative_number = -2.1;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Sign(positive_number, +1.1), positive_number);
	ASSERT_DOUBLE_EQ(Sign(positive_number, -1.1), -1.0 * positive_number);
	ASSERT_DOUBLE_EQ(Sign(negative_number, +1.1), -1.0 * negative_number);
	ASSERT_DOUBLE_EQ(Sign(negative_number, -1.1), negative_number);
}

TEST(TestSpecialFunctions, TestStepFunction)
{
	// ARRANGE
	double positive_number = 10.3;
	double zero			   = 0.0;
	double negative_number = -2.1;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(StepFunction(positive_number), 1.0);
	ASSERT_DOUBLE_EQ(StepFunction(negative_number), 0.0);
	ASSERT_DOUBLE_EQ(StepFunction(zero), 1.0);
}

TEST(TestSpecialFunctions, TestRound)
{
	// ARRANGE
	double pi = 3.14159265359;
	double no = 5.6563782e-17;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Round(pi, 1), 3.0);
	ASSERT_DOUBLE_EQ(Round(pi, 2), 3.1);
	ASSERT_DOUBLE_EQ(Round(pi, 3), 3.14);
	ASSERT_DOUBLE_EQ(Round(pi, 4), 3.142);
	ASSERT_DOUBLE_EQ(Round(pi, 5), 3.1416);
	ASSERT_DOUBLE_EQ(Round(pi, 6), 3.14159);
	ASSERT_DOUBLE_EQ(Round(pi, 7), 3.141593);
	ASSERT_DOUBLE_EQ(Round(no, 1), 6.0e-17);
	ASSERT_DOUBLE_EQ(Round(no, 2), 5.7e-17);
	ASSERT_DOUBLE_EQ(Round(no, 3), 5.66e-17);
	ASSERT_DOUBLE_EQ(Round(no, 4), 5.656e-17);
	ASSERT_DOUBLE_EQ(Round(no, 5), 5.6564e-17);
	ASSERT_DOUBLE_EQ(Round(no, 6), 5.65638e-17);
	ASSERT_DOUBLE_EQ(Round(no, 7), 5.656378e-17);
}

TEST(TestSpecialFunctions, TestRoundVector)
{
	// ARRANGE
	double pi  = 3.14159265359;
	double no  = 5.6563782e-17;
	double neg = -19.1457;
	Vector vec({pi, no, 2.0 * pi, neg});
	// ACT & ASSERT
	for(unsigned int i = 0; i < vec.Size(); i++)
		ASSERT_DOUBLE_EQ(Round(vec, 1)[i], Vector({3.0, 6.0e-17, 6.0, -20.0})[i]);
	for(unsigned int i = 0; i < vec.Size(); i++)
		ASSERT_DOUBLE_EQ(Round(vec, 2)[i], Vector({3.1, 5.7e-17, 6.3, -19.0})[i]);
	for(unsigned int i = 0; i < vec.Size(); i++)
		ASSERT_DOUBLE_EQ(Round(vec, 3)[i], Vector({3.14, 5.66e-17, 6.28, -19.1})[i]);
}

TEST(TestSpecialFunctions, TestRoundMatrix)
{
	// ARRANGE
	double pi  = 3.14159265359;
	double no  = 5.6563782e-17;
	double neg = -19.1457;
	Matrix M({{pi, no}, {neg, 2.0 * pi}});
	// ACT & ASSERT
	for(unsigned int i = 0; i < M.Rows(); i++)
		for(unsigned int j = 0; j < M.Columns(); j++)
			ASSERT_DOUBLE_EQ(Round(M, 1)[i][j], Matrix({{3.0, 6.0e-17}, {-20.0, 6.0}})[i][j]);
	for(unsigned int i = 0; i < M.Rows(); i++)
		for(unsigned int j = 0; j < M.Columns(); j++)
			ASSERT_DOUBLE_EQ(Round(M, 2)[i][j], Matrix({{3.1, 5.7e-17}, {-19.0, 6.3}})[i][j]);
	for(unsigned int i = 0; i < M.Rows(); i++)
		for(unsigned int j = 0; j < M.Columns(); j++)
			ASSERT_DOUBLE_EQ(Round(M, 3)[i][j], Matrix({{3.14, 5.66e-17}, {-19.1, 6.28}})[i][j]);
}

TEST(TestSpecialFunctions, TestRelativeDifference)
{
	// ARRANGE
	double number_1			   = 1.3;
	double number_2			   = 1.0;
	double relative_difference = 3.0 / 13.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Relative_Difference(number_1, number_2), relative_difference);
	ASSERT_DOUBLE_EQ(Relative_Difference(number_1, number_1), 0.0);
}

TEST(TestSpecialFunctions, TestFloats_Equal)
{
	// ARRANGE
	double float_1 = 1.0000000000;
	double float_2 = 1.0000000002;
	// ACT & ASSERT
	ASSERT_TRUE(Floats_Equal(float_1, float_2, 1e-6));
	ASSERT_FALSE(Floats_Equal(float_1, float_2, 1e-10));
}

// 2. Special functions
// 2.1 Gamma functions
TEST(TestSpecialFunctions, TestFactorial)
{
	// ARRANGE
	std::vector<int> numbers			= {0, 1, 2, 3, 4, 5, 10, 20, 170};
	std::vector<double> correct_results = {1, 1, 2, 6, 24, 120, 3628800, 2432902008176640000, 7.25741561530799896739672821113e306};
	// ACT & ASSERT
	for(unsigned int i = 0; i < numbers.size(); i++)
		ASSERT_DOUBLE_EQ(Factorial(numbers[i]), correct_results[i]);
}

TEST(TestSpecialFunctions, TestBinomialCoefficient)
{
	// ACT & ASSERT
	ASSERT_EQ(Binomial_Coefficient(3, 5), 0);
	ASSERT_EQ(Binomial_Coefficient(5, 0), 1);
	ASSERT_EQ(Binomial_Coefficient(5, 1), 5);
	ASSERT_EQ(Binomial_Coefficient(5, 2), 10);
	ASSERT_EQ(Binomial_Coefficient(5, 3), 10);
	ASSERT_EQ(Binomial_Coefficient(5, 4), 5);
	ASSERT_EQ(Binomial_Coefficient(5, 5), 1);
	ASSERT_EQ(Binomial_Coefficient(50, 10), 10272278170);
}

TEST(TestSpecialFunctions, TestGammaLn)
{
	// ARRANGE
	double tolerance = 1.0e-10;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaLn(1.0), 0.0);
	ASSERT_NEAR(GammaLn(3.0), log(2.0), tolerance);
	ASSERT_NEAR(GammaLn(3.0 / 2.0), log(sqrt(M_PI) / 2.0), tolerance);
	ASSERT_NEAR(GammaLn(4.0), log(6.0), tolerance);
	ASSERT_NEAR(GammaLn(8.0), log(5040.0), tolerance);
}

TEST(TestSpecialFunctions, TestGamma)
{
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Gamma(2.0), 1.0);
	ASSERT_DOUBLE_EQ(Gamma(7.5), 135135 * sqrt(M_PI) / 128.0);
	for(unsigned int k = 0; k < 11; k++)
		ASSERT_NEAR(Gamma(k + 1), Factorial(k), 1.0e-6);
}

TEST(TestSpecialFunctions, TestIncompleteGamma)
{
	// ARRANGE
	double s = 3.5;
	// ACT & ASSERT
	ASSERT_NEAR(Upper_Incomplete_Gamma(1.0, 4.0), 16.0 / exp(1.0), 1.0e-12);
	for(int x = 0; x < 10; x++)
		ASSERT_NEAR(Upper_Incomplete_Gamma(x, s) + Lower_Incomplete_Gamma(x, s), Gamma(s), 1.0e-12);
}

TEST(TestSpecialFunctions, TestGammaQ)
{
	// ARRANGE
	double e = exp(1.0);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaQ(0.0, 3.0), 1.0);
	ASSERT_DOUBLE_EQ(GammaQ(1.0, 1.0), 1.0 / e);
	ASSERT_DOUBLE_EQ(GammaQ(1.0, 3.0), 5.0 / 2.0 / e);
	ASSERT_NEAR(GammaQ(5.0, 2.0), 6.0 / e / e / e / e / e, 1.0e-4);
	ASSERT_NEAR(GammaQ(110.0, 110.0), 0.48732, 1.0e-6);
}

TEST(TestSpecialFunctions, TestGammaP)
{
	// ARRANGE
	double e = exp(1.0);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaP(1.0, 3.0), 1.0 - 5.0 / (2.0 * e));
}

TEST(TestSpecialFunctions, TestInvGammaP)
{
	// ARRANGE
	double p = 0.7;
	double a = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaP(Inv_GammaP(p, a), a), p);
}

TEST(TestSpecialFunctions, TestInvGammaQ)
{
	// ARRANGE
	double q = 0.7;
	double a = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaQ(Inv_GammaQ(q, a), a), q);
}

// 2.2 Scalar spherical Harmonics

TEST(TestSpecialFunctions, TestSphericalHarmonics)
{
	// ARRANGE
	int l						= 5;
	int m						= 3;
	double theta				= 25 * deg;
	double phi					= 13 * deg;
	std::complex<double> result = -0.129726 - 1.0i * 0.10505;
	double tol					= 1.0e-6;
	// ACT & ASSERT
	EXPECT_NEAR(Spherical_Harmonics(1, 0, theta, phi).real(), 0.5 * sqrt(3.0 / M_PI) * cos(theta), tol);

	EXPECT_NEAR(Spherical_Harmonics(l, m, theta, phi).real(), result.real(), tol);
	EXPECT_NEAR(Spherical_Harmonics(l, m, theta, phi).imag(), result.imag(), tol);

	EXPECT_NEAR(Spherical_Harmonics(l, m, 0.0, 0.0).imag(), 0.0, tol);
	EXPECT_NEAR(Spherical_Harmonics(l, m, 0.0, 0.0).imag(), 0.0, tol);

	EXPECT_NEAR(Spherical_Harmonics(l, 2 * l, theta, phi).imag(), 0.0, tol);
	EXPECT_NEAR(Spherical_Harmonics(l, 2 * l, theta, phi).imag(), 0.0, tol);
}

// 2.3 Scalar spherical Harmonics
TEST(TestSpecialFunctions, TestVectorialSphericalHarmonicsY)
{
	// ARRANGE
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Y	= Vector_Spherical_Harmonics_Y(l, m, theta, phi);
	auto Ym = Vector_Spherical_Harmonics_Y(l, -m, theta, phi);
	// ASSERT
	for(int i = 0; i < 3; i++)
	{
		EXPECT_NEAR(Ym[i].real(), std::pow(-1.0, m) * std::conj(Y[i]).real(), 1.0e-6);
		EXPECT_NEAR(Ym[i].imag(), std::pow(-1.0, m) * std::conj(Y[i]).imag(), 1.0e-6);
	}
}

TEST(TestSpecialFunctions, TestVectorialSphericalHarmonicsPsi)
{
	// ARRANGE
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Psi  = Vector_Spherical_Harmonics_Psi(l, m, theta, phi);
	auto Psim = Vector_Spherical_Harmonics_Psi(l, -m, theta, phi);
	// ASSERT
	for(int i = 0; i < 3; i++)
	{
		EXPECT_NEAR(Psim[i].real(), std::pow(-1.0, m) * std::conj(Psi[i]).real(), 1.0e-6);
		EXPECT_NEAR(Psim[i].imag(), std::pow(-1.0, m) * std::conj(Psi[i]).imag(), 1.0e-6);
	}
}

TEST(TestSpecialFunctions, TestVectorialSphericalHarmonicsOrthogonality)
{
	// ARRANGE
	int l		 = 5;
	int m		 = 3;
	double theta = 25 * deg;
	double phi	 = 13 * deg;
	// ACT
	auto Y	 = Vector_Spherical_Harmonics_Y(l, m, theta, phi);
	auto Psi = Vector_Spherical_Harmonics_Psi(l, m, theta, phi);
	// ASSERT
	EXPECT_NEAR((Y[0] * Psi[0] + Y[1] * Psi[1] + Y[2] * Psi[2]).real(), 0.0, 1e-16);
	EXPECT_NEAR((Y[0] * Psi[0] + Y[1] * Psi[1] + Y[2] * Psi[2]).imag(), 0.0, 1e-16);
}

TEST(TestSpecialFunctions, TestVectorialSphericalHarmonicsYFailure)
{
	// ACT & ASSERT
	ASSERT_DEATH({ auto comp = VSH_Y_Component(5, 5, 3, 5, 3); }, "");
}

TEST(TestSpecialFunctions, TestVectorialSphericalHarmonicsPsiFailure)
{
	// ACT & ASSERT
	ASSERT_DEATH({ auto comp = VSH_Psi_Component(5, 5, 3, 5, 3); }, "");
}

// 2.4 Other special functions
TEST(TestSpecialFunctions, TestDawsonIntegral)
{
	// ARRANGE
	std::vector<double> xs		= {0.1, 0.5, 2.3, 10.7, 22.0};
	std::vector<double> results = {0.099336, 0.424436, 0.249053, 0.0469358, 0.0227508};
	// ACT & ASSERT
	for(unsigned int i = 0; i < xs.size(); i++)
		EXPECT_NEAR(Dawson_Integral(xs[i]), results[i], 1e-6);
}

TEST(TestSpecialFunctions, TestErfi)
{
	// ARRANGE
	std::vector<double> xs		= {0.01, 0.2, 4.2, 8.4, -2.0 / 3.0, 0.15, 0.45, 2.35, 10.24, 22.9};
	std::vector<double> results = {0.0112842, 0.228721, 6.34555e6, 2.97919e29, -0.880276, 0.170535, 0.544232, 68.3326, 1.91557e44, 1.38157e226};
	double tol					= 1.0e-5;
	// ACT & ASSERT
	for(unsigned int i = 0; i < xs.size(); i++)
		EXPECT_NEAR(Erfi(xs[i]), results[i], std::fabs(tol * results[i]));
}

TEST(TestSpecialFunctions, TestInvErf)
{
	// ARRANGE
	double y = 0.45;
	// ACT & ASSERT
	ASSERT_NEAR(erf(Inv_Erf(y)), y, 1.0e-6);
}