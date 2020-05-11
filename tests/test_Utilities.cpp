#include "gtest/gtest.h"

#include "Utilities.hpp"

//1. Progress bar

//2. Import and export data from files
	TEST(TestUtilities, TestExportImportList)
	{
		// ARRANGE
		double unit = 2.3;
		std::vector<double> list = {1.0 * unit, 2.0 * unit, 3.0 * unit, 4.0 * unit, 5.0 * unit};
		// ACT
		Export_List("test_list.txt", list, unit);
		std::vector<double> imported_list = Import_List("test_list.txt", unit);
		// ASSERT
		ASSERT_EQ( imported_list.size(), list.size() );
		for(unsigned int i = 0; i < list.size(); i++)
			ASSERT_DOUBLE_EQ( imported_list[i], list[i] );
	}

	TEST(TestUtilities, TestExportImportTable)
	{
		// ARRANGE
		std::vector<double> units = {2.3, 3.1};
		std::vector<std::vector<double>> table = {{11.0 * units[0], 12.0 * units[1]},{21.0 * units[0], 22.0 * units[1]}};
		// ACT
		Export_Table("test_table.txt", table, units);
		std::vector<std::vector<double>> imported_table = Import_Table("test_table.txt",units);
		// ASSERT
		ASSERT_EQ( imported_table.size(), table.size() );
		for(unsigned int i = 0; i < table.size(); i++)
		{
			ASSERT_EQ( imported_table[i].size(), table[i].size() );
			for(unsigned int j = 0; j < table[i].size(); j++)
				ASSERT_DOUBLE_EQ( imported_table[i][j], table[i][j]);			
		}		
	}

	double func(double x) {return x*x;}
	TEST(TestUtilities, TestExportImportFunction)
	{
		// ARRANGE
		double xMin = 1.0;
		double xMax = 3.0;
		unsigned int steps = 3;
		std::vector<std::vector<double>> correct_values = {{1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}};
		// ACT
		Export_Function("test_function.txt", func, xMin, xMax, steps);
		std::vector<std::vector<double>> imported_values = Import_Table("test_function.txt");
		// ASSERT
		ASSERT_EQ( imported_values.size(), correct_values.size() );
		for(unsigned int i = 0; i < correct_values.size(); i++)
		{
			ASSERT_EQ( imported_values[i].size(), correct_values[i].size() );
			for(unsigned int j = 0; j < correct_values[i].size(); j++)
				ASSERT_DOUBLE_EQ( imported_values[i][j], correct_values[i][j]);			
		}		
	}
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

	TEST(TestUtilities, TestLogSpace)
		{
		// ARRANGE
		double min = 1.0e1;
		double max = 1.0e10;
		unsigned int steps = 10;
		std::vector<double> list_correct = {1.0e1,1.0e2,1.0e3,1.0e4,1.0e5,1.0e6,1.0e7,1.0e8,1.0e9,1.0e10};
		// ACT
		std::vector<double> list = Log_Space(min,max,steps);
		// ASSERT
		ASSERT_EQ( list.size(), list_correct.size() );
		for(unsigned int i = 0; i < list_correct.size(); i++)
			ASSERT_NEAR( list[i], list_correct[i] , 1.0e-3);
	}