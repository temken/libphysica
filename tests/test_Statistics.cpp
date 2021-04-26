#include "gtest/gtest.h"

#include "libphysica/Statistics.hpp"

using namespace libphysica;

//1. Distributions
//1.1. Uniform distribution
TEST(TestStatistics, TestPDFUniform)
{
	// ARRANGE
	double x_Min	   = 10.0;
	double x_Max	   = 12.0;
	double x		   = 11.3;
	double PDF_correct = 1.0 / (x_Max - x_Min);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Uniform(x, x_Min, x_Max), PDF_correct);
}

TEST(TestStatistics, TestCDFUniform)
{
	// ARRANGE
	double x_Min = 10.0;
	double x_Max = 12.0;
	double x	 = x_Min + (x_Max - x_Min) / 3.0;
	double x_Mid = (x_Max + x_Min) / 2.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Uniform(x_Min, x_Min, x_Max), 0.0);
	ASSERT_NEAR(CDF_Uniform(x, x_Min, x_Max), 1.0 / 3.0, 1.0e-10);
	ASSERT_DOUBLE_EQ(CDF_Uniform(x_Mid, x_Min, x_Max), 0.5);
	ASSERT_DOUBLE_EQ(CDF_Uniform(x_Max, x_Min, x_Max), 1.0);
}

//1.2 Normal distribution
TEST(TestStatistics, TestPDFGauss)
{
	// ARRANGE
	double mu		 = 2.5;
	double sigma	 = 1.2;
	double tolerance = 1.0e-8;
	auto f			 = std::bind(PDF_Gauss, std::placeholders::_1, mu, sigma);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Gauss(mu, mu, sigma), 1.0 / sqrt(2.0 * M_PI) / sigma);
	ASSERT_NEAR(Integrate(f, mu - sigma, mu + sigma, tolerance), 0.682689492137, tolerance);
	ASSERT_NEAR(Integrate(f, mu - 2.0 * sigma, mu + 2.0 * sigma, tolerance), 0.954499736104, tolerance);
	ASSERT_NEAR(Integrate(f, mu - 3.0 * sigma, mu + 3.0 * sigma, tolerance), 0.997300203937, tolerance);
}

TEST(TestStatistics, TestCDFGauss)
{
	// ARRANGE
	double mu		   = 2.5;
	double sigma	   = 1.2;
	double p_one_sigma = 0.682689492137;
	// ACT & ASSERT
	ASSERT_NEAR(CDF_Gauss(mu - sigma, mu, sigma), 0.5 - p_one_sigma / 2.0, 1.0e-10);
	ASSERT_DOUBLE_EQ(CDF_Gauss(mu, mu, sigma), 0.5);
	ASSERT_NEAR(CDF_Gauss(mu + sigma, mu, sigma), 0.5 + p_one_sigma / 2.0, 1.0e-10);
}

TEST(TestStatistics, TestQuantileGauss)
{
	// ARRANGE
	double mu	 = 2.5;
	double sigma = 1.2;
	double p	 = 0.75;
	// ACT & ASSERT
	ASSERT_NEAR(CDF_Gauss(Quantile_Gauss(p, mu, sigma), mu, sigma), p, 1.0e-8);
}

//1.3 Binomial distribution
TEST(TestStatistics, TestPMFBinomial)
{
	// ARRANGE
	unsigned int trials = 10;
	double p			= 0.25;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PMF_Binomial(trials, p, 3), 0.25028228759765625);
	ASSERT_DOUBLE_EQ(PMF_Binomial(trials, p, 6), 0.0162220001220703125);
}

TEST(TestStatistics, TestCDFBinomial)
{
	// ARRANGE
	unsigned int trials = 10;
	double p			= 0.25;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Binomial(trials, p, 3), 0.775875091552734375);
	ASSERT_DOUBLE_EQ(CDF_Binomial(trials, p, 6), 0.996494293212890625);
}

//1.4 Poission distribution
TEST(TestStatistics, TestPMFPoisson)
{
	// ARRANGE
	double expectation_value = M_PI;
	// ACT & ASSERT
	ASSERT_NEAR(PMF_Poisson(expectation_value, 3), 0.2233171177, 1.0e-10);
	ASSERT_NEAR(PMF_Poisson(expectation_value, 5), 0.1102025804, 1.0e-10);
}

TEST(TestStatistics, TestCDFPoisson)
{
	// ARRANGE
	double expectation_value = M_PI;
	// ACT & ASSERT
	ASSERT_NEAR(CDF_Poisson(expectation_value, 3), 0.6155437030, 1.0e-10);
	ASSERT_NEAR(CDF_Poisson(expectation_value, 5), 0.9011391375, 1.0e-10);
}

TEST(TestStatistics, TestInvCDFPoisson)
{
	// ARRANGE
	unsigned int x = 0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Inv_CDF_Poisson(x, 0.1), log(10));
	ASSERT_DOUBLE_EQ(Inv_CDF_Poisson(x, 0.05), log(20));
	ASSERT_DOUBLE_EQ(Inv_CDF_Poisson(x, 0.01), log(100));
}
//1.5 Chi-square distribution
TEST(TestStatistics, TestPDFChiSquare)
{
	// ARRANGE
	double dof = 2.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Chi_Square(-1.0, dof), 0.0);
	ASSERT_DOUBLE_EQ(PDF_Chi_Square(0.0, dof), 0.0);
	ASSERT_DOUBLE_EQ(PDF_Chi_Square(1.0, 0.0), 0.0);
	ASSERT_DOUBLE_EQ(PDF_Chi_Square(0.4, dof), 1.0 / 2.0 / exp(0.2));
}

TEST(TestStatistics, TestCDFChiSquare)
{
	// ARRANGE
	double dof = 2.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Chi_Square(-1.0, dof), 0.0);
	ASSERT_DOUBLE_EQ(CDF_Chi_Square(0.0, dof), 0.0);
	ASSERT_DOUBLE_EQ(CDF_Chi_Square(0.0, 0.0), 1.0);
	ASSERT_DOUBLE_EQ(CDF_Chi_Square(1.0, 0.0), 1.0);
	ASSERT_DOUBLE_EQ(CDF_Chi_Square(0.4, dof), 1.0 - 1.0 / exp(0.2));
}

TEST(TestStatistics, TestPDFChiBarSquare)
{
	// ARRANGE
	std::vector<double> weights = {0.6, 0.3, 0.07, 0.03};
	double x					= 1.5;
	double correct_pdf			= weights[1] * PDF_Chi_Square(x, 1.0) + weights[2] * PDF_Chi_Square(x, 2.0) + weights[3] * PDF_Chi_Square(x, 3.0);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Chi_Bar_Square(x, {0.0, 1.0}), PDF_Chi_Square(x, 1.0));
	ASSERT_DOUBLE_EQ(PDF_Chi_Bar_Square(x, weights), correct_pdf);
}

TEST(TestStatistics, TestCDFChiBarSquare)
{
	// ARRANGE
	std::vector<double> weights = {0.6, 0.3, 0.07, 0.03};
	double x					= 1.5;
	double correct_cdf			= weights[0] + weights[1] * CDF_Chi_Square(x, 1.0) + weights[2] * CDF_Chi_Square(x, 2.0) + weights[3] * CDF_Chi_Square(x, 3.0);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Chi_Bar_Square(x, {0.0, 1.0}), CDF_Chi_Square(x, 1.0));
	ASSERT_DOUBLE_EQ(CDF_Chi_Bar_Square(x, weights), correct_cdf);
	ASSERT_DOUBLE_EQ(CDF_Chi_Bar_Square(0.0, weights), weights[0]);
	ASSERT_DOUBLE_EQ(CDF_Chi_Bar_Square(-1.0, weights), 0.0);
}

//1.6 Exponential distribution
TEST(TestStatistics, TestPDFExponential)
{
	// ARRANGE
	double mean = 1.5;
	double x	= 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Exponential(x, mean), exp(-3.0 / mean) / mean);
	ASSERT_DOUBLE_EQ(PDF_Exponential(0, mean), 1.0 / mean);
	ASSERT_DOUBLE_EQ(PDF_Exponential(-1.0, mean), 0.0);
}

TEST(TestStatistics, TestCDFExponential)
{
	// ARRANGE
	double mean = 1.5;
	double x	= 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Exponential(x, mean), 1.0 - exp(-3.0 / mean));
	ASSERT_DOUBLE_EQ(CDF_Exponential(0.0, mean), 0.0);
	ASSERT_DOUBLE_EQ(CDF_Exponential(-1.0, mean), 0.0);
}

//1.6 Exponential distribution
TEST(TestStatistics, TestPDFMaxwellBoltzmann)
{
	// ARRANGE
	double a = 1.5;
	double x = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(PDF_Maxwell_Boltzmann(x, a), 8.0 * sqrt(2.0 / M_PI) / 3.0 * exp(-2.0));
	ASSERT_DOUBLE_EQ(PDF_Maxwell_Boltzmann(-1.0, a), 0.0);
	ASSERT_DOUBLE_EQ(PDF_Maxwell_Boltzmann(0.0, a), 0.0);
}

TEST(TestStatistics, TestCDFMaxwellBoltzmann)
{
	// ARRANGE
	double a = 1.5;
	double x = 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(CDF_Maxwell_Boltzmann(x, a), erf(sqrt(2)) - 2.0 * sqrt(2.0 / M_PI) * exp(-2.0));
	ASSERT_DOUBLE_EQ(CDF_Maxwell_Boltzmann(-1.0, a), 0.0);
	ASSERT_DOUBLE_EQ(CDF_Maxwell_Boltzmann(0.0, a), 0.0);
	ASSERT_DOUBLE_EQ(CDF_Maxwell_Boltzmann(100. * a, a), 1.0);
}

//2. Likelihoods
// extern double Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background = 0.0);
TEST(TestStatistics, TestLikelihoodPoisson)
{
	// ARRANGE
	double mu	   = 3.0;
	unsigned int x = 4;
	double B	   = 0.2;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Likelihood_Poisson(mu, x, B), 8192.0 / 1875.0 / exp(16.0 / 5.0));
}
// extern double Log_Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background = 0.0);
TEST(TestStatistics, TestLogLikelihoodPoisson)
{
	// ARRANGE
	double mu	   = 3.0;
	unsigned int x = 4;
	double B	   = 0.2;
	double mu_tot  = mu + B;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Log_Likelihood_Poisson(mu, x, B), -mu_tot + 4 * log(mu_tot) - log(24));
}
// extern double Likelihood_Poisson_Binned(std::vector<double> N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned = {});
TEST(TestStatistics, TestLikelihoodPoissonBinned)
{
	// ARRANGE
	std::vector<double> mu				= {3.3, 2.1, 0.5};
	std::vector<unsigned long int> data = {0, 3, 2};
	double llh_correct					= 1.0;
	for(unsigned int i = 0; i < mu.size(); i++)
		llh_correct *= Likelihood_Poisson(mu[i], data[i]);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Likelihood_Poisson_Binned(mu, data), llh_correct);
}
// extern double Log_Likelihood_Poisson_Binned(std::vector<double> N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned = {});
TEST(TestStatistics, TestLogLikelihoodPoissonBinned)
{
	// ARRANGE
	std::vector<double> mu				= {3.3, 2.1, 0.5};
	std::vector<unsigned long int> data = {0, 3, 2};
	double log_llh_correct				= 0.0;
	for(unsigned int i = 0; i < mu.size(); i++)
		log_llh_correct += Log_Likelihood_Poisson(mu[i], data[i]);
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Log_Likelihood_Poisson_Binned(mu, data), log_llh_correct);
}
//3. Sampling random numbers
// 3.1 Sample from specific distribution
TEST(TestStatistics, SampleUniform)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	double xMin = -3.0;
	double xMax = 4.0;
	// ACT & ASSERT
	for(int i = 0; i < 10; i++)
	{
		double x = Sample_Uniform(PRNG, xMin, xMax);
		ASSERT_GT(x, xMin);
		ASSERT_LT(x, xMax);
	}
}

// extern int Sample_Poisson(std::mt19937& PRNG, double expectation_value);
TEST(TestStatistics, TestSamplePoisson1)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	double mu = 4.0;
	// ACT & ASSERT
	for(int i = 0; i < 10; i++)
		ASSERT_GE(Sample_Poisson(PRNG, mu), 0.0);
}

TEST(TestStatistics, TestSamplePoisson2)
{
	// ARRANGE
	unsigned int fixed_seed = 15;
	std::mt19937 PRNG(fixed_seed);
	std::vector<unsigned int> result = {10, 14, 12};

	std::vector<double> expectation_values = {10.3, 12.5, 8.7};
	// ACT
	std::vector<unsigned int> samples = Sample_Poisson(PRNG, expectation_values);
	// ASSERT
	ASSERT_EQ(samples.size(), expectation_values.size());
	for(unsigned int i = 0; i < samples.size(); i++)
		EXPECT_EQ(samples[i], result[i]);
}

//3.2 General sampling algorithms
// extern double Rejection_Sampling(const std::function<double(double)>& PDF,double xMin,double xMax,double yMax,std::mt19937& PRNG);
TEST(TestStatistics, TestRejectionSampling)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	double lambda = 3.0;
	double xMin	  = 0.0;
	double xMax	  = 20.0;
	auto pdf	  = std::bind(PDF_Exponential, std::placeholders::_1, lambda);
	// ACT & ASSERT
	for(int i = 0; i < 10; i++)
	{
		double x = Rejection_Sampling(pdf, xMin, xMax, pdf(xMin), PRNG);
		ASSERT_GT(x, xMin);
		ASSERT_LT(x, xMax);
	}
}

// extern double Inverse_Transform_Sampling(const std::function<double(double)>& cdf,double xMin,double xMax,std::mt19937& PRNG);
TEST(TestStatistics, TestInverseTransformSampling)
{
	// ARRANGE
	std::random_device rd;
	std::mt19937 PRNG(rd());
	double mu	 = 0.0;
	double sigma = 2.0;
	double xMin	 = -100.0;
	double xMax	 = 100.0;
	auto cdf	 = std::bind(CDF_Gauss, std::placeholders::_1, mu, sigma);
	// ACT & ASSERT
	for(int i = 0; i < 10; i++)
	{
		double x = Inverse_Transform_Sampling(cdf, xMin, xMax, PRNG);
		ASSERT_GT(x, xMin);
		ASSERT_LT(x, xMax);
	}
}
//4. Data point with statistical weight
TEST(TestStatistics, TestDataPoint)
{
	// ARRANGE
	DataPoint P1(1.0);
	DataPoint P2(2.0);
	DataPoint P3(1.0);
	// ACT & ASSERT
	ASSERT_LT(P1, P2);
	ASSERT_GT(P2, P1);
	ASSERT_EQ(P1, P3);
}

//5. Basic data analysis
TEST(TestStatistics, TestBasicDataAnalysis)
{
	// ARRANGE
	std::vector<double> data	 = {0.5, 1.0, 1.2, 1.5};
	std::vector<DataPoint> data2 = {DataPoint(1.0), DataPoint(2.0)};
	// ACT & ASSERT
	EXPECT_DOUBLE_EQ(Arithmetic_Mean(data), 21.0 / 20.0);
	EXPECT_DOUBLE_EQ(Median(data), 11.0 / 10.0);
	EXPECT_DOUBLE_EQ(Variance(data), 53.0 / 300.0);
	EXPECT_DOUBLE_EQ(Standard_Deviation(data), sqrt(53.0 / 300.0));
	EXPECT_DOUBLE_EQ(Weighted_Average(data2)[0], 1.5);
}

//6. Kernel density estimation
// extern Interpolation Perform_KDE(std::vector<DataPoint> data,double xMin,double xMax,double bw = 0);
TEST(TestStatistics, TestPerformKDE)
{
	// ARRANGE
	std::vector<DataPoint> data = {DataPoint(13.0, 1.0), DataPoint(9.0, 0.5), DataPoint(11.0, 1.0), DataPoint(16.0, 1.0)};
	double xMin					= 0.0;
	double xMax					= 30.0;
	double tolerance			= 1.0e-6;
	//ACT
	Interpolation kde = Perform_KDE(data, xMin, xMax);
	double norm		  = Integrate(kde, xMin, 29.99999, tolerance);
	//ASSERT
	ASSERT_GT(kde(10.0), 0.0);
	ASSERT_NEAR(norm, 1.0, tolerance);
}
