#include "gtest/gtest.h"

#include "Statistics.hpp"

//1. Distributions
	TEST(TestStatistics, TestPDFUniform)
	{
		// ARRANGE
		double x_Min = 10.0;
		double x_Max = 12.0;
		double x = 11.3;
		double PDF_correct = 1.0 / (x_Max - x_Min);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( PDF_Uniform(x, x_Min, x_Max), PDF_correct );
	}

	TEST(TestStatistics, TestCDFUniform)
	{
		// ARRANGE
		double x_Min = 10.0;
		double x_Max = 12.0;
		double x = x_Min + (x_Max - x_Min) / 3.0;
		double x_Mid = (x_Max + x_Min) / 2.0;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( CDF_Uniform(x_Min, x_Min, x_Max), 0.0 );
		ASSERT_NEAR( CDF_Uniform(x, x_Min, x_Max), 1.0/3.0 ,1.0e-10);
		ASSERT_DOUBLE_EQ( CDF_Uniform(x_Mid, x_Min, x_Max), 0.5 );
		ASSERT_DOUBLE_EQ( CDF_Uniform(x_Max, x_Min, x_Max), 1.0 );
	}

	TEST(TestStatistics, TestPDFGauss)
	{
		// ARRANGE
		double mu = 2.5;
		double sigma = 1.2;
		double tolerance = 1.0e-8;
		auto f = std::bind(PDF_Gauss, std::placeholders::_1, mu, sigma);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( PDF_Gauss(mu, mu, sigma), 1.0/sqrt(2.0*M_PI)/sigma);
		ASSERT_NEAR( Integrate(f, mu-sigma, mu+sigma, tolerance), 0.682689492137, tolerance);
		ASSERT_NEAR( Integrate(f, mu-2.0*sigma, mu+2.0*sigma, tolerance), 0.954499736104, tolerance);
		ASSERT_NEAR( Integrate(f, mu-3.0*sigma, mu+3.0*sigma, tolerance), 0.997300203937, tolerance);
	}

	TEST(TestStatistics, TestCDFGauss)
	{
		// ARRANGE
		double mu = 2.5;
		double sigma = 1.2;
		double p_one_sigma = 0.682689492137;
		// ACT & ASSERT
		ASSERT_NEAR( CDF_Gauss(mu - sigma, mu, sigma), 0.5 - p_one_sigma/2.0, 1.0e-10);
		ASSERT_DOUBLE_EQ( CDF_Gauss(mu, mu, sigma), 0.5);
		ASSERT_NEAR( CDF_Gauss(mu + sigma, mu, sigma), 0.5 + p_one_sigma/2.0, 1.0e-10);
	}

	TEST(TestStatistics, TestQuantileGauss)
	{
		// ARRANGE
		double mu = 2.5;
		double sigma = 1.2;		
		double p = 0.75;
		// ACT & ASSERT
		ASSERT_NEAR( CDF_Gauss(Quantile_Gauss(p, mu, sigma), mu, sigma), p, 1.0e-8);
	}


//2. Likelihoods

//3. Sampling random numbers

//4. Data point with statistical weight
	TEST(TestStatistics, TestDataPoint)
	{
		// ARRANGE
		DataPoint P1(1.0);
		DataPoint P2(2.0);
		DataPoint P3(1.0);
		// ACT & ASSERT
		ASSERT_LT( P1 , P2 );
		ASSERT_GT( P2 , P1 );
		ASSERT_EQ( P1 , P3 );
	}

//5. Basic data analysis
	TEST(TestStatistics, TestBasicDataAnalysis)
	{
		// ARRANGE
		std::vector<double> data = {0.5, 1.0, 1.2, 1.5};
		std::vector<DataPoint> data2 = {DataPoint(1.0), DataPoint(2.0)};
		// ACT & ASSERT
		EXPECT_DOUBLE_EQ( Arithmetic_Mean(data), 21.0/20.0);
		EXPECT_DOUBLE_EQ( Median(data), 11.0/10.0);
		EXPECT_DOUBLE_EQ( Variance(data), 53.0/300.0);
		EXPECT_DOUBLE_EQ( Standard_Deviation(data), sqrt(53.0/300.0));
		EXPECT_DOUBLE_EQ( Weighted_Average(data2)[0], 1.5);

	}

//6. Kernel density estimation
