#include "gtest/gtest.h"

#include "Linear_Algebra.hpp"

TEST(TestVector, TestSize)
{
	//ARRANGE
	unsigned int dim = 4;
	Vector v(dim);
	//ACT & ASSERT
	ASSERT_EQ( v.Size() ,  dim);
}

TEST(TestVector, TestDot)
{
	//ARRANGE
	Vector v1({1,2,3});
	Vector v2(3, 2.0);
	double correct_result = 12.0;
	//ACT & ASSERT
	ASSERT_DOUBLE_EQ(v1.Dot(v2) ,  correct_result);
}