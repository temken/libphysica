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

// 1.2 1D integration with boost functions
TEST(TestNumerics, TestIntegrateBoost)
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

TEST(TestNumerics, TestIntegrate2D)
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

TEST(TestNumerics, TestIntegrate3D)
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

TEST(TestNumerics, TestIntegrate3DVector)
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

TEST(TestNumerics, TestIntegrateMCBruteForce)
{
	// ARRANGE
	double tolerance		   = 1.0e-3;
	std::vector<double> region = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
	double result			   = std::pow((exp(1.0) - 1.0) / exp(1), 3.0);
	auto func				   = [](const std::vector<double>& x, const double wgt) {
		 return exp(-x[0] - x[1] - x[2]);
	};
	// ACT & ASSERT
	EXPECT_NEAR(Integrate_MC(func, region, 1000000, "Brute force"), result, tolerance);
}

TEST(TestNumerics, TestIntegrateMCVegas)
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