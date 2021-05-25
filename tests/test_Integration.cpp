#include "gtest/gtest.h"

#include <cmath>
#include <functional>
#include <random>

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