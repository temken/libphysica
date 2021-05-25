#ifndef __Integration_hpp__
#define __Integration_hpp__

#include <functional>
#include <string>
#include <vector>

namespace libphysica
{

// 1. Integration
// 1.1 One-dimensional integration via adaptive Simpson method
extern double Find_Epsilon(std::function<double(double)> func, double a, double b, double precision);
extern double Integrate(std::function<double(double)> func, double a, double b, double epsilon, int maxRecursionDepth = 20);

// 2. Multidimensional MC integration
extern double Integrate_MC(std::function<double(std::vector<double>&)> fxn, std::vector<double>& region, double& standard_deviation, unsigned int function_evaluations = 10000, const std::string& method = "Vegas");

}	// namespace libphysica

#endif