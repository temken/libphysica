#include "gtest/gtest.h"

#include "Statistics.hpp"

//1. Distributions
	TEST(TestDistributions, TestPDFUniform)
	{
		// ARRANGE
		double x_Min = 10.0;
		double x_Max = 12.0;
		double x = 11.3;
		double PDF_correct = 1.0 / (x_Max - x_Min);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( PDF_Uniform(x,x_Min, x_Max), PDF_correct );
	}
//2. Likelihoods

//3. Sampling random numbers

//4. Data point with statistical weight

//5. Basic data analysis

//6. Kernel density estimation
