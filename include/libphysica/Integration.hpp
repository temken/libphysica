#ifndef __Integration_hpp__
#define __Integration_hpp__

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "libphysica/Linear_Algebra.hpp"

namespace libphysica
{

// 1. Integration
// 1.1 One-dimensional integration via adaptive Simpson method
extern double Find_Epsilon(std::function<double(double)> func, double a, double b, double precision);
extern double Integrate(std::function<double(double)> func, double a, double b, double epsilon, int maxRecursionDepth = 20);

// 1.2 Gauss Legendre quadrature
extern std::vector<std::vector<double>> Compute_Gauss_Legendre_Roots_and_Weights(unsigned int n, double x_min = -1.0, double x_max = 1.0);
extern double Integrate_Gauss_Legendre(std::function<double(double)> func, double a, double b, unsigned int sample_points = 30);
extern double Integrate_Gauss_Legendre(std::function<double(double)> func, std::vector<std::vector<double>> roots_and_weights);
extern double Integrate_Gauss_Legendre(std::vector<double> function_values, std::vector<std::vector<double>> roots_and_weights);

// 1.3 1D integration with boost functions
extern double Integrate(std::function<double(double)> func, double a, double b, const std::string& method = "Gauss-Legendre", int method_parameter = 0);

// 2. Multidimensional integration
// 2.1 Multidimensional integration via nesting 1D integration
extern double Integrate_2D(std::function<double(double, double)> func, double x1, double x2, double y1, double y2, const std::string& method = "Gauss-Legendre", int method_parameter = 0);
extern double Integrate_3D(std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2, const std::string& method = "Gauss-Legendre", int method_parameter = 0);
extern double Integrate_3D(std::function<double(Vector)> func, double r1, double r2, double costheta_1 = -1.0, double costheta_2 = 1.0, double phi_1 = 0.0, double phi_2 = 2.0 * M_PI, const std::string& method = "Gauss-Legendre", int method_parameter = 0);

// 2.2 Monte Carlo Integration
extern double Integrate_MC(std::function<double(std::vector<double>&, const double)> func, std::vector<double>& region, const int ncalls = 10000, const std::string& method = "Vegas");

}	// namespace libphysica

#endif