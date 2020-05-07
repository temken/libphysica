#include "gtest/gtest.h"

#include "Utilities.hpp"

//1. Progress bar

//2. Import and export data from files

//3. Create list with equi-distant numbers in log-space
	TEST(TestUtilities, TestLinearSpace)
	{
		// ARRANGE
		double min = 1.0;
		double max = 10.0;
		unsigned int steps = 10;
		std::vector<double> list_correct = {1,2,3,4,5,6,7,8,9,10};
		// ACT
		std::vector<double> list = Linear_Space(min,max,steps);
		// ASSERT
		ASSERT_EQ( list.size(), list_correct.size() );
		for(unsigned int i = 0; i < list_correct.size(); i++)
			ASSERT_DOUBLE_EQ( list[i], list_correct[i] );
	}