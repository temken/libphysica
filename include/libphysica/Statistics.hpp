#ifndef __Statistics_hpp_
#define __Statistics_hpp_

#include <random>
#include <vector>

#include "Numerics.hpp"

namespace libphysica
{

//1. Distributions
//1.1. Uniform distribution
extern double PDF_Uniform(double x, double x_min, double x_max);
extern double CDF_Uniform(double x, double x_min, double x_max);

//1.2 Normal distribution
extern double PDF_Gauss(double x, double mu, double sigma);
extern double CDF_Gauss(double x, double mu, double sigma);
extern double Quantile_Gauss(double p, double mu, double sigma);

//1.3 Binomial distribution
extern double PMF_Binomial(unsigned int trials, double p, unsigned int x);
extern double CDF_Binomial(unsigned int trials, double p, unsigned int x);

//1.4 Poission distribution
extern double PMF_Poisson(double expected_events, unsigned int events);
extern double CDF_Poisson(double expectation_value, unsigned int observed_events);
extern double Inv_CDF_Poisson(unsigned int observed_events, double cdf);   //Solves the CDF = cdf for mu

//1.5 Chi-square distribution
extern double PDF_Chi_Square(double x, double dof);
extern double CDF_Chi_Square(double x, double dof);
extern double PDF_Chi_Bar_Square(double x, std::vector<double> weights);
extern double CDF_Chi_Bar_Square(double x, std::vector<double> weights);

//1.6 Exponential distribution
extern double PDF_Exponential(double x, double mean);
extern double CDF_Exponential(double x, double mean);

//1.7 Maxwell-Boltzmann distribution
extern double PDF_Maxwell_Boltzmann(double x, double a);
extern double CDF_Maxwell_Boltzmann(double x, double a);

//2. Likelihoods
extern double Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background = 0.0);
extern double Log_Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background = 0.0);
extern double Likelihood_Poisson_Binned(const std::vector<double>& N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned = {});
extern double Log_Likelihood_Poisson_Binned(const std::vector<double>& N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned = {});

//3. Sampling random numbers
//3.1 Sample from specific distribution
extern double Sample_Uniform(std::mt19937& PRNG, double x_min = 0.0, double x_max = 1.0);
extern unsigned int Sample_Poisson(std::mt19937& PRNG, double expectation_value);
extern std::vector<unsigned int> Sample_Poisson(std::mt19937& PRNG, const std::vector<double>& expectation_values);

//3.2 General sampling algorithms
extern double Rejection_Sampling(const std::function<double(double)>& PDF, double xMin, double xMax, double yMax, std::mt19937& PRNG);
extern double Inverse_Transform_Sampling(const std::function<double(double)>& cdf, double xMin, double xMax, std::mt19937& PRNG);

//4. Data point with statistical weight
struct DataPoint
{
	double value;
	double weight;
	//Constructors:
	DataPoint(double v = 0.0, double w = 1.0);
};
bool operator<(const DataPoint& lhs, const DataPoint& rhs);
bool operator>(const DataPoint& lhs, const DataPoint& rhs);
bool operator==(const DataPoint& lhs, const DataPoint& rhs);
std::ostream& operator<<(std::ostream& output, const DataPoint& dp);

//5. Basic data analysis
extern double Arithmetic_Mean(const std::vector<double>& data);
extern double Median(std::vector<double>& data);
extern double Variance(const std::vector<double>& data);
extern double Standard_Deviation(const std::vector<double>& data);
extern std::vector<double> Weighted_Average(std::vector<DataPoint>& data);

//6. Kernel density estimation
extern Interpolation Perform_KDE(std::vector<DataPoint> data, double xMin, double xMax, double bw = 0);

}	// namespace libphysica

#endif
