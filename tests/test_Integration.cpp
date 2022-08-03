#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <random>
#include <string>
#include <vector>

#include "libphysica/Integration.hpp"
// #include "libphysica/Statistics.hpp"
// #include "libphysica/Utilities.hpp"

using namespace libphysica;

// 1. Integration
// 1.1 One-dimensional integration via adaptive Simpson method
TEST(TestIntegration, TestIntegrate)
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

// 1.2 Gauss Legendre quadrature
TEST(TestIntegration, TestGaussLegendreRoots)
{
	// ARRANGE
	int n = 4;
	std::vector<double> roots {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
	std::vector<double> weights {0.34785485, 0.65214515, 0.65214515, 0.34785485};

	// ACT
	auto roots_and_weights = Compute_Gauss_Legendre_Roots_and_Weights(n);

	// ASSERT
	for(int i = 0; i < n; i++)
	{
		EXPECT_FLOAT_EQ(roots_and_weights[i][0], roots[i]);
		EXPECT_FLOAT_EQ(roots_and_weights[i][1], weights[i]);
	}
}

TEST(TestIntegration, TestGaussLegendre)
{
	// ARRANGE
	std::function<double(double)> func = [](double x) {
		return x * x * x - 2.0 * x * x + x - 1.0;
	};
	double a = 3.4, b = 4.5;
	unsigned int n = 4;
	// ACT
	double integral = Integrate_Gauss_Legendre(func, a, b, n);
	// ASSERT
	EXPECT_NEAR(integral, 37.8049, 1.0e-4);
}

// 1.3 1D integration with boost functions
TEST(TestIntegration, TestIntegrateBoost)
{
	// ARRANGE
	double tolerance				   = 1.0e-8;
	std::function<double(double)> func = [](double x) {
		return 1.0 / x;
	};
	std::function<double(double)> func2 = [](double x) {
		return sin(x);
	};
	std::vector<std::string> methods = {"Trapezoidal", "Gauss-Legendre", "Gauss-Kronrod"};
	// ACT & ASSERT
	for(auto& method : methods)
	{
		EXPECT_NEAR(Integrate(func, 2.0, 3.0, method), log(3.0 / 2.0), tolerance);
		EXPECT_NEAR(Integrate(func2, -1.0, 1.0, method), 0.0, tolerance);
	}
}

TEST(TestIntegration, TestIntegrationLimits)
{
	// ARRANGE
	double tolerance				   = 1.0e-8;
	std::function<double(double)> func = [](double x) {
		return 1.0 / x;
	};
	double a						 = 2.0;
	double b						 = 3.0;
	std::vector<std::string> methods = {"Trapezoidal", "Gauss-Legendre", "Gauss-Kronrod"};

	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(Integrate(func, a, b, tolerance), -1.0 * Integrate(func, b, a, tolerance));
	for(auto& method : methods)
		EXPECT_DOUBLE_EQ(Integrate(func, a, b, method), -1.0 * Integrate(func, b, a, method));
}

TEST(TestIntegration, TestIntegrate2DNested)
{
	// ARRANGE
	double tolerance = 1.0e-6;
	auto func		 = [](double x, double y) {
		   return exp(-x * x - y * y);
	};
	std::vector<std::string> methods = {"Trapezoidal", "Gauss-Legendre", "Gauss-Kronrod"};
	// ACT & ASSERT
	for(auto& method : methods)
		EXPECT_NEAR(Integrate_2D(func, -1.0, 1.0, -1.0, 1.0, method), M_PI * erf(1.0) * erf(1.0), tolerance);
}

TEST(TestIntegration, TestIntegrate2DMC)
{
	// ARRANGE
	double tolerance_MC	   = 1.0e-2;
	double tolerance_Vegas = 1.0e-4;
	auto func			   = [](double x, double y) {
		 return exp(-x * x - y * y);
	};
	// ACT & ASSERT
	EXPECT_NEAR(Integrate_2D(func, -1.0, 1.0, -1.0, 1.0, "Monte-Carlo"), M_PI * erf(1.0) * erf(1.0), tolerance_MC);
	EXPECT_NEAR(Integrate_2D(func, -1.0, 1.0, -1.0, 1.0, "Vegas"), M_PI * erf(1.0) * erf(1.0), tolerance_Vegas);
}

TEST(TestIntegration, TestIntegrate3DNested)
{
	// ARRANGE
	double tolerance = 1.0e-6;
	auto func		 = [](double x, double y, double z) {
		   return x + y + z;
	};
	std::vector<std::string> methods = {"Trapezoidal", "Gauss-Legendre", "Gauss-Kronrod"};
	// ACT & ASSERT
	for(auto& method : methods)
		EXPECT_NEAR(Integrate_3D(func, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, method), 1.5, tolerance);
}

TEST(TestIntegration, TestIntegrate3DMC)
{
	// ARRANGE
	double tolerance_MC	   = 1.0e-2;
	double tolerance_Vegas = 1.0e-4;
	auto func			   = [](double x, double y, double z) {
		 return x + y + z;
	};
	// ACT & ASSERT
	EXPECT_NEAR(Integrate_3D(func, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, "Monte-Carlo"), 1.5, tolerance_MC);
	EXPECT_NEAR(Integrate_3D(func, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, "Vegas"), 1.5, tolerance_Vegas);
}

TEST(TestIntegration, TestIntegrate3DVector)
{
	// ARRANGE
	double tolerance = 1.0e-6;
	double result	 = 625.0 * M_PI;
	auto func		 = [](Vector vector) {
		   return vector.Norm();
	};
	std::vector<std::string> methods = {"Trapezoidal", "Gauss-Legendre", "Gauss-Kronrod"};
	// ACT & ASSERT
	for(auto& method : methods)
		EXPECT_NEAR(Integrate_3D(func, 0.0, 5.0, -1.0, 1.0, 0.0, 2.0 * M_PI, method), result, tolerance * result);
}

TEST(TestIntegration, TestIntegrateMCBruteForce)
{
	// ARRANGE
	double tolerance		   = 1.0e-3;
	std::vector<double> region = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double result			   = std::pow((exp(1.0) - 1.0) / exp(1), 3.0);
	auto func				   = [](const std::vector<double>& x, const double wgt) {
		 return exp(-x[0] - x[1] - x[2]);
	};
	// ACT & ASSERT
	EXPECT_NEAR(Integrate_MC(func, region, 1000000, "Monte-Carlo"), result, tolerance);
}

TEST(TestIntegration, TestIntegrateMCVegas)
{
	// ARRANGE
	double tolerance		   = 1.0e-4;
	std::vector<double> region = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double result			   = std::pow((exp(1.0) - 1.0) / exp(1), 3.0);
	auto func				   = [](const std::vector<double>& x, const double wgt) {
		 return exp(-x[0] - x[1] - x[2]);
	};
	// ACT & ASSERT
	EXPECT_NEAR(Integrate_MC(func, region, 1000000, "Vegas"), result, tolerance);
}