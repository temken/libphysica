#include "gtest/gtest.h"

#include "Numerics.hpp"


//1. Simple functions
	TEST(TestNumerics, TestSign)
	{
		// ARRANGE
		double positive_number = 10.3;
		double negative_number = -2.1;
		// ACT & ASSERT
		ASSERT_EQ(Sign(positive_number), 1);
		ASSERT_EQ(Sign(negative_number), -1);
	}
//2. Special functions

//2.1 Gamma functions

//2.2 Other special functions

//3. Integration

//3.1 One-dimensional integration via adaptive Simpson method 

//4. Interpolation

//4.1 One-dimensional interpolation

//5. Root finding

//6. Minimization

//6.1 Multi-dimensional
