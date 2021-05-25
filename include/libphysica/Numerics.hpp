#ifndef __Numerics_hpp_
#define __Numerics_hpp_

#include <functional>
#include <string>
#include <vector>

namespace libphysica
{

// 1. Interpolation
// 1.1 One-dimensional interpolation (Steffen splines)
class Interpolation
{
  private:
	unsigned int N;
	std::vector<double> x_values;
	std::vector<double> function_values;
	double prefactor;

	//Steffen coefficients
	std::vector<double> a, b, c, d;
	void Compute_Steffen_Coefficients();

	//Locate j in array such that x[j]<x<x[j+1].
	unsigned int jLast;
	bool correlated_calls;	 // if successive calls are correlated, then the hunt method can be faster.
	unsigned int Bisection(double x, int jLeft, int jRight);
	unsigned int Hunt(double x);

  public:
	std::vector<double> domain;

	unsigned int Locate(double x);

	//Constructors
	Interpolation();
	explicit Interpolation(const std::vector<double>& arg_values, const std::vector<double>& func_values, double x_dim = -1.0, double f_dim = -1.0);
	explicit Interpolation(const std::vector<std::vector<double>>& data, double x_dim = -1.0, double f_dim = -1.0);

	void Set_Prefactor(double factor);
	void Multiply(double factor);

	double Interpolate(double x);
	double Derivative(double x, unsigned int deriv = 1);
	double Integrate(double x_1, double x_2);
	double operator()(double x)
	{
		return Interpolate(x);
	};

	void Save_Function(std::string filename, unsigned int points);
};

// 1.2 Two-dimensional interpolation (bilinear interpolation)
class Interpolation_2D
{
  private:
	unsigned int N_x, N_y;
	std::vector<double> x_values;
	std::vector<double> y_values;
	std::vector<std::vector<double>> function_values;
	double prefactor;
	Interpolation x_int, y_int;	  //Dummy 1D interpolations, necessary for the Locate() functionality.

  public:
	std::vector<std::vector<double>> domain;
	Interpolation_2D();
	Interpolation_2D(std::vector<double> x_val, std::vector<double> y_val, std::vector<std::vector<double>> func_values, double x_dim = -1.0, double y_dim = -1.0, double f_dim = -1.0);
	Interpolation_2D(std::vector<std::vector<double>> data_table, double x_dim = -1.0, double y_dim = -1.0, double f_dim = -1.0);

	double Interpolate(double x, double y);

	double operator()(double x, double y)
	{
		return Interpolate(x, y);
	};

	void Set_Prefactor(double factor);
	void Multiply(double factor);

	void Save_Function(std::string filename, unsigned int x_points, unsigned int y_points = 0);
};

// 2. Root finding
extern double Find_Root(std::function<double(double)> func, double xLeft, double xRight, double epsilon);

// 3. Minimization
// 3.1 Multi-dimensional
struct Minimization
{
	const double ftol;
	int nfunc;	 //The number of function evaluations.
	int mpts;
	int ndim;
	double fmin;										//Function value at the minimum.
	std::vector<double> y;								//Function values at the vertices of the simplex.
	std::vector<std::vector<double>> current_simplex;	// p
	explicit Minimization(const double ftoll)
	: ftol(ftoll) {}   //The constructor argument ftoll is the fractional convergence tolerance to be achieved in the function value (n.b.!).

	//Multidimensional minimization of the function or functor func(x), where x[0..ndim-1] is a vector in ndim dimensions, by the downhill simplex method of Nelder and Mead. The initial simplex is specified as in equation (10.5.1) by a point[0..ndim-1] and a constant displacement del along each coordinate direction. Returned is the location of the minimum.
	std::vector<double> minimize(std::vector<double>& starting_point, const double delta, std::function<double(std::vector<double>)> func);

	//Alternative interface that takes different displacements dels[0..ndim-1] in different di- rections for the initial simplex.
	std::vector<double> minimize(std::vector<double>& starting_point, std::vector<double>& deltas, std::function<double(std::vector<double>)> func);

	//Most general interface: initial simplex specified by the matrix pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors that are the vertices of the starting simplex.
	std::vector<double> minimize(std::vector<std::vector<double>>& pp, std::function<double(std::vector<double>)> func);

	void get_psum(std::vector<std::vector<double>>& p, std::vector<double>& psum);	 //Utility function.

	double amotry(std::vector<std::vector<double>>& p, std::vector<double>& y, std::vector<double>& psum, const int ihi, const double fac, std::function<double(std::vector<double>)> func);
};

}	// namespace libphysica

#endif