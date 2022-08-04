#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <random>

#include "libphysica/Numerics.hpp"
#include "libphysica/Statistics.hpp"
#include "libphysica/Utilities.hpp"

using namespace libphysica;

// 1. Interpolation
// 1.1 One-dimensional interpolation
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

TEST(TestNumerics, TestInterpolation1DDomainExtrapolation)
{
	// ARRANGE
	std::vector<std::vector<double>> data = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}};
	// ACT
	Interpolation interpolation(data);
	// ASSERT
	ASSERT_GT(interpolation(3.001), 0.0);
	ASSERT_LT(interpolation(-0.005), 0.0);
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

TEST(TestNumerics, TestInterpolation1DIntegrate)
{
	// ARRANGE
	std::vector<double> x_values = Linear_Space(-5, 5, 100);
	std::vector<double> y_values;
	for(auto& x : x_values)
		y_values.push_back(x * x * x);
	Interpolation interpol(x_values, y_values);
	double tol = 1.0e-5;
	// ACT & ASSERT
	ASSERT_NEAR(interpol.Integrate(-4.4, 4.4), 0.0, tol);
	ASSERT_NEAR(interpol.Integrate(1, 3.5), 2385.0 / 64, tol);
	ASSERT_NEAR(interpol.Integrate(3.5, 1), -2385.0 / 64, tol);
	ASSERT_NEAR(interpol.Integrate(0, M_PI), M_PI * M_PI * M_PI * M_PI / 4.0, tol);
}

TEST(TestNumerics, TestInterpolation1DLocalMinimum)
{
	// ARRANGE
	std::vector<double> x_values = Linear_Space(-5, 5, 10000);
	std::vector<double> y_values, y_values_2;
	for(auto& x : x_values)
		y_values.push_back(x * x * x);
	Interpolation interpol(x_values, y_values);
	for(auto& x : x_values)
		y_values_2.push_back(x * x);
	Interpolation interpol_2(x_values, y_values_2);
	double tol = 1e-6;
	// ACT & ASSERT
	EXPECT_NEAR(interpol.Local_Minimum(-1.0, 1.0), -1.0, tol);
	EXPECT_NEAR(interpol.Local_Minimum(2.0, 2.3), 8.0, tol);
	EXPECT_NEAR(interpol_2.Local_Minimum(-2.0, 2.0), 0.0, tol);
	ASSERT_DEATH(interpol.Local_Minimum(2.4, 2.1), "");
}

TEST(TestNumerics, TestInterpolation1DLocalMaximum)
{
	// ARRANGE
	std::vector<double> x_values = Linear_Space(-5, 5, 10000);
	std::vector<double> y_values, y_values_2;
	for(auto& x : x_values)
		y_values.push_back(x * x * x);
	Interpolation interpol(x_values, y_values);
	for(auto& x : x_values)
		y_values_2.push_back(sin(x));
	Interpolation interpol_2(x_values, y_values_2);
	double tol = 1e-6;
	// ACT & ASSERT
	EXPECT_NEAR(interpol.Local_Maximum(-1.0, 1.0), 1.0, tol);
	EXPECT_NEAR(interpol.Local_Maximum(2.0, 2.3), 2.3 * 2.3 * 2.3, tol);
	EXPECT_NEAR(interpol_2.Local_Maximum(0.0, 1.0), sin(1), tol);
	EXPECT_NEAR(interpol_2.Local_Maximum(0.0, 3.0), 1.0, tol);
}

TEST(TestNumerics, TestInterpolation1DGlobalMinimum)
{
	// ARRANGE
	std::vector<double> x_values = Linear_Space(-5, 5, 10000);
	std::vector<double> y_values;
	for(auto& x : x_values)
		y_values.push_back(x * x - 5.0);
	Interpolation interpol(x_values, y_values);
	double tol = 1e-6;
	// ACT & ASSERT
	ASSERT_NEAR(interpol.Global_Minimum(), -5.0, tol);
}

TEST(TestNumerics, TestInterpolation1DGlobalMaximum)
{
	// ARRANGE
	std::vector<double> x_values = Linear_Space(-5, 5, 10000);
	std::vector<double> y_values;
	for(auto& x : x_values)
		y_values.push_back(-std::pow(x - 2.3, 2) + M_PI);
	Interpolation interpol(x_values, y_values);
	double tol = 1e-6;
	// ACT & ASSERT
	ASSERT_NEAR(interpol.Global_Maximum(), M_PI, tol);
}

// 1.2 Two-dimensional interpolation (bilinear interpolation)
TEST(TestNumerics, TestInterpolation2dDefaultConstructor)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 100;
	// ACT
	Interpolation_2D interpolation;
	// ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, -1.0, 1.0);
		double y = Sample_Uniform(PRNG, -1.0, 1.0);
		ASSERT_DOUBLE_EQ(interpolation(x, y), 0.0);
	}
}

TEST(TestNumerics, TestInterpolation2dListConstructor)
{
	// ARRANGE
	std::vector<double> x_list				 = {-5.0, 0.0, 5.0};
	std::vector<double> y_list				 = {-5.0, 0.0, 5.0};
	std::vector<std::vector<double>> f_table = {{-5.0, 0.0, 5.0}, {-5.0, 0.0, 5.0}, {-5.0, 0.0, 5.0}};

	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 100;
	double tolerance	= 1.0e-8;
	// ACT
	Interpolation_2D interpolation(x_list, y_list, f_table);
	// ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, interpolation.domain[0][0], interpolation.domain[0][1]);
		double y = Sample_Uniform(PRNG, interpolation.domain[1][0], interpolation.domain[1][1]);
		double f = y;
		ASSERT_NEAR(interpolation(x, y), f, tolerance);
	}
}

TEST(TestNumerics, TestInterpolation2dTableConstructor)
{
	// ARRANGE
	std::vector<std::vector<double>> data_table = {
		{-5.0, -5.0, -5.0},
		{-5.0, 0.0, 0.0},
		{-5.0, 5.0, 5.0},
		{0.0, -5.0, -5.0},
		{0.0, 0.0, 0.0},
		{0.0, 5.0, 5.0},
		{5.0, -5.0, -5.0},
		{5.0, 0.0, 0.0},
		{5.0, 5.0, 5.0}};

	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 100;
	double tolerance	= 1.0e-8;
	// ACT
	Interpolation_2D interpolation(data_table);
	// ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, interpolation.domain[0][0], interpolation.domain[0][1]);
		double y = Sample_Uniform(PRNG, interpolation.domain[1][0], interpolation.domain[1][1]);
		double f = y;
		ASSERT_NEAR(interpolation(x, y), f, tolerance);
	}
}

TEST(TestNumerics, TestInterpolation2dGlobalMinimum)
{
	// ARRANGE
	std::vector<std::vector<double>> data_table = {
		{-5.0, -5.0, -53.0},
		{-5.0, 0.0, 0.0},
		{-5.0, 5.0, 5.0},
		{0.0, -5.0, -5.0},
		{0.0, 0.0, 0.0},
		{0.0, 5.0, 5.0},
		{5.0, -5.0, -5.0},
		{5.0, 0.0, 0.0},
		{5.0, 5.0, 8.0}};
	Interpolation_2D interpolation(data_table);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(interpolation.Global_Minimum(), -53.0);
}

TEST(TestNumerics, TestInterpolation2dGlobalMinimum2)
{
	// ARRANGE
	std::vector<std::vector<double>> data_table = {
		{-5.0, -5.0, -53.0},
		{-5.0, 0.0, 0.0},
		{-5.0, 5.0, 5.0},
		{0.0, -5.0, -5.0},
		{0.0, 0.0, 0.0},
		{0.0, 5.0, 5.0},
		{5.0, -5.0, -5.0},
		{5.0, 0.0, 0.0},
		{5.0, 5.0, 8.0}};
	Interpolation_2D interpolation(data_table);

	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 10000;

	// ACT & ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, interpolation.domain[0][0], interpolation.domain[0][1]);
		double y = Sample_Uniform(PRNG, interpolation.domain[1][0], interpolation.domain[1][1]);
		EXPECT_TRUE(interpolation(x, y) > interpolation.Global_Minimum());
	}
}

TEST(TestNumerics, TestInterpolation2dGlobalMaximum)
{
	// ARRANGE
	std::vector<double> x_list				 = {-5.0, 0.0, 5.0};
	std::vector<double> y_list				 = {-5.0, 0.0, 5.0};
	std::vector<std::vector<double>> f_table = {{-5.0, 0.0, 5.0}, {-5.0, 13.6, 5.0}, {-5.0, 0.0, 5.0}};
	Interpolation_2D interpolation(x_list, y_list, f_table);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(interpolation.Global_Maximum(), 13.6);
}

TEST(TestNumerics, TestInterpolation2dGlobalMaximum2)
{
	// ARRANGE
	std::vector<double> x_list				 = {-5.0, 0.0, 5.0};
	std::vector<double> y_list				 = {-5.0, 0.0, 5.0};
	std::vector<std::vector<double>> f_table = {{-5.0, 0.0, 5.0}, {-5.0, 13.6, 5.0}, {-5.0, 0.0, 5.0}};
	Interpolation_2D interpolation(x_list, y_list, f_table);

	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 10000;

	// ACT & ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, interpolation.domain[0][0], interpolation.domain[0][1]);
		double y = Sample_Uniform(PRNG, interpolation.domain[1][0], interpolation.domain[1][1]);
		EXPECT_TRUE(interpolation(x, y) < interpolation.Global_Maximum());
	}
}

TEST(TestNumerics, TestInterpolation2dUnits)
{
	// ARRANGE
	std::vector<double> x_list				 = {-5.0, 0.0, 5.0};
	std::vector<double> y_list				 = {-5.0, 0.0, 5.0};
	std::vector<std::vector<double>> f_table = {{-5.0, 0.0, 5.0}, {-5.0, 0.0, 5.0}, {-5.0, 0.0, 5.0}};
	double dim_x							 = 3;
	double dim_y							 = 5;
	double dim_f							 = 7;
	std::random_device rd;
	std::mt19937 PRNG(rd());
	unsigned int trials = 100;
	double tolerance	= 1.0e-8;
	// ACT
	Interpolation_2D interpolation(x_list, y_list, f_table, dim_x, dim_y, dim_f);
	// ASSERT
	for(unsigned int i = 0; i < trials; i++)
	{
		double x = Sample_Uniform(PRNG, -5.0, 5.0);
		double y = Sample_Uniform(PRNG, -5.0, 5.0);
		double f = y;
		ASSERT_NEAR(interpolation(x * dim_x, y * dim_y), f * dim_f, tolerance);
	}
}

// 2. Root finding
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

// 3. Minimization
// 3.1 One-dimensional
double minimize_function_1d_1(double x)
{
	double sig = std::pow(2.2, -1.0 / 6);
	return std::pow(sig / x, 12) - std::pow(sig / x, 6);
}
double minimize_function_1d_2(double x)
{
	return cos(x);
}
TEST(TestNumerics, TestMinimization1D)
{
	// ARRANGE
	double minimum_1 = 0.9842404697785967;
	// ACT & ASSERT
	EXPECT_NEAR(Find_Minimum(minimize_function_1d_1, 0.0, 2.0), minimum_1, 1.0e-8);
	EXPECT_NEAR(Find_Minimum(minimize_function_1d_2, 2.0, 4.0), M_PI, 1.0e-8);
}

double maximize_function_1d(double x)
{
	double C = 1.9;
	return -std::pow(x - C, 2.0);
}

TEST(TestNumerics, TestMaximization1D)
{
	// ARRANGE
	double maximum = 1.9;
	// ACT & ASSERT
	EXPECT_NEAR(Find_Maximum(maximize_function_1d, 0.0, 10.), maximum, 1.0e-8);
}

// 3.2 Multi-dimensional
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