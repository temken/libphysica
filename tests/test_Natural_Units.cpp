#include "gtest/gtest.h"

#include "Natural_Units.hpp"

// 1. SI-prefixes

//2. Units

//3. Physical constants

//4. Astronomical parameters

//5. Transform quantities to a given dimension
	TEST(TestInUnits, TestInUnitsScalar)
	{
		//ARRANGE
		double c = 1.0;
		double speed_of_light_mps = 299792458.0;
		//ACT & ASSERT
		ASSERT_DOUBLE_EQ( In_Units(c, meter/sec), speed_of_light_mps );
	}
	
//6. Simple physics functions

