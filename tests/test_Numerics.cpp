#include "gtest/gtest.h"

#include <cmath>
#include <functional>

#include "Numerics.hpp"

using namespace libphysica;

//1. Simple functions
TEST(TestNumerics, TestSign)
{
	// ARRANGE
	double positive_number = 10.3;
	double negative_number = -2.1;
	// ACT & ASSERT
	ASSERT_EQ(Sign(positive_number), 1);
	ASSERT_EQ(Sign(negative_number), -1);
}

TEST(TestNumerics, TestSign2)
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

TEST(TestNumerics, TestStepFunction)
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

TEST(TestNumerics, TestRound)
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

TEST(TestNumerics, TestRelativeDifference)
{
	// ARRANGE
	double number_1			   = 1.3;
	double number_2			   = 1.0;
	double relative_difference = 3.0 / 13.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Relative_Difference(number_1, number_2), relative_difference);
	ASSERT_DOUBLE_EQ(Relative_Difference(number_1, number_1), 0.0);
}

TEST(TestNumerics, TestFloats_Equal)
{
	// ARRANGE
	double float_1 = 1.0000000000;
	double float_2 = 1.0000000002;
	// ACT & ASSERT
	ASSERT_TRUE(Floats_Equal(float_1, float_2, 1e-6));
	ASSERT_FALSE(Floats_Equal(float_1, float_2, 1e-10));
}

//2. Special functions
//2.1 Gamma functions
TEST(TestNumerics, TestFactorial)
{
	// ARRANGE
	std::vector<int> numbers			= {0, 1, 2, 3, 4, 5, 10, 20, 170};
	std::vector<double> correct_results = {1, 1, 2, 6, 24, 120, 3628800, 2432902008176640000, 7.25741561530799896739672821113e306};
	// ACT & ASSERT
	for(unsigned int i = 0; i < numbers.size(); i++)
		ASSERT_DOUBLE_EQ(Factorial(numbers[i]), correct_results[i]);
}

TEST(TestNumerics, TestBinomialCoefficient)
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

TEST(TestNumerics, TestGammaLn)
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

TEST(TestNumerics, TestGamma)
{
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Gamma(2.0), 1.0);
	ASSERT_DOUBLE_EQ(Gamma(7.5), 135135 * sqrt(M_PI) / 128.0);
	for(unsigned int k = 0; k < 11; k++)
		ASSERT_NEAR(Gamma(k + 1), Factorial(k), 1.0e-6);
}

TEST(TestNumerics, TestIncompleteGamma)
{
	// ARRANGE
	double s = 3.5;
	// ACT & ASSERT
	ASSERT_NEAR(Upper_Incomplete_Gamma(1.0, 4.0), 16.0 / exp(1.0), 1.0e-12);
	for(int x = 0; x < 10; x++)
		ASSERT_NEAR(Upper_Incomplete_Gamma(x, s) + Lower_Incomplete_Gamma(x, s), Gamma(s), 1.0e-12);
}

TEST(TestNumerics, TestGammaQ)
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

TEST(TestNumerics, TestGammaP)
{
	// ARRANGE
	double e = exp(1.0);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaP(1.0, 3.0), 1.0 - 5.0 / (2.0 * e));
}

TEST(TestNumerics, TestInvGammaP)
{
	// ARRANGE
	double p = 0.7;
	double a = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaP(Inv_GammaP(p, a), a), p);
}

TEST(TestNumerics, TestInvGammaQ)
{
	// ARRANGE
	double q = 0.7;
	double a = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(GammaQ(Inv_GammaQ(q, a), a), q);
}

//2.2 Other special functions
TEST(TestNumerics, TestInvErf)
{
	// ARRANGE
	double y = 0.45;
	// ACT & ASSERT
	ASSERT_NEAR(erf(Inv_Erf(y)), y, 1.0e-6);
}

//3. Integration
//3.1 One-dimensional integration via adaptive Simpson method
TEST(TestNumerics, TestIntegrate)
{
	// ARRANGE
	double tolerance				   = 1.0e-8;
	std::function<double(double)> func = [](double x) {
		return 1.0 / x;
	};
	std::function<double(double)> func2 = [](double x) {
		return sin(x);
	};
	// ACT & ASSERT
	ASSERT_NEAR(Integrate(func, 2.0, 3.0, tolerance), log(3.0 / 2.0), tolerance);
	ASSERT_NEAR(Integrate(func2, -1.0, 1.0, tolerance), 0.0, tolerance);
}

//4. Interpolation
//4.1 One-dimensional interpolation
TEST(TestNumerics, TestInterpolation1D1)
{
	// ARRANGE
	std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
	std::vector<double> f = {0.0, 1.0, 2.0, 3.0};
	// ACT
	Interpolation interpolation(x, f);
	// ASSERT
	for(unsigned int i = 0; i < x.size(); i++)
		ASSERT_DOUBLE_EQ(interpolation(x[i]), f[i]);
	for(unsigned int i = 0; i < x.size() - 1; i++)
		ASSERT_DOUBLE_EQ(interpolation(x[i] + 0.5), f[i] + 0.5);
}

TEST(TestNumerics, TestInterpolation1D2)
{
	// ARRANGE
	std::vector<std::vector<double>> data = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}};
	// ACT
	Interpolation interpolation(data);
	// ASSERT
	for(unsigned int i = 0; i < data.size(); i++)
		ASSERT_DOUBLE_EQ(interpolation(data[i][0]), data[i][1]);
	for(unsigned int i = 0; i < data.size() - 1; i++)
		ASSERT_DOUBLE_EQ(interpolation(data[i][0] + 0.5), data[i][1] + 0.5);
}

TEST(TestNumerics, TestInterpolation1DDerivative)
{
	// ARRANGE
	double offset						  = 2.0;
	std::vector<std::vector<double>> data = {{-3.0, 9.0 + offset}, {-2.0, 4.0 + offset}, {-1.0, 1.0 + offset}, {0.0, offset}, {1.0, 1.0 + offset}, {2.0, 4.0 + offset}, {3.0, 9.0 + offset}};
	double x							  = 1.5;
	// ACT
	Interpolation interpolation(data);
	// ASSERT
	EXPECT_DOUBLE_EQ(interpolation(x), x * x + offset);
	EXPECT_DOUBLE_EQ(interpolation.Derivative(x, 0), x * x + offset);
	EXPECT_DOUBLE_EQ(interpolation.Derivative(x), 2.0 * x);
	EXPECT_DOUBLE_EQ(interpolation.Derivative(x, 2), 2.0);
	EXPECT_DOUBLE_EQ(interpolation.Derivative(x, 3), 0.0);
	EXPECT_DOUBLE_EQ(interpolation.Derivative(x, 4), 0.0);
}
//5. Root finding
double find_root_func(double x)
{
	return x * x - 2.0;
}
TEST(TestNumerics, TestFindRoot)
{
	// ARRANGE
	double epsilon = 1.0e-6;
	// ACT & ASSERT
	ASSERT_NEAR(Find_Root(find_root_func, 0.0, 2.0, epsilon), sqrt(2.0), epsilon);
	ASSERT_NEAR(Find_Root(find_root_func, -2.0, 0.0, epsilon), -sqrt(2.0), epsilon);
}

//6. Minimization
//6.1 Multi-dimensional
double minimize_function(std::vector<double> args)
{
	return (args[0] - M_PI) * (args[0] - M_PI) + (args[1] - exp(1.0)) * (args[1] - exp(1.0));
}
TEST(TestNumerics, TestMinimization)
{
	// ARRANGE
	double tolerance = 1.0e-6;
	Minimization min(tolerance);
	std::vector<double> starting_point = {0.0, 0.0};
	double delta					   = 0.5;
	// ACT
	std::vector<double> minimum = min.minimize(starting_point, delta, minimize_function);
	// ASSERT
	ASSERT_NEAR(minimum[0], M_PI, tolerance);
	ASSERT_NEAR(minimum[1], exp(1.0), tolerance);
}