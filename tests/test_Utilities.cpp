#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "libphysica/Utilities.hpp"

using namespace libphysica;

// 1. Progress bar
TEST(TestUtilities, TestTimeDisplay)
{
	// ARRANGE
	// ACT
	std::vector<std::string> results = {"[00m:00s:001ms]", "[00m:00s:010ms]", "[00m:00s:100ms]", "[00m:01s:000ms]", "[00m:10s:000ms]", "[01m:40s:000ms]", "[16m:40s:000ms]", "[02h:46m:40s]", "[01d:03h:46m]", "[01w:04d:13h]", "[16w:03d:17h]", "[03y:08w:05d]", "[31y:35w:06d]", "[316y:45w:06d]"};
	// ACT & ASSERT
	for(int i = -3; i < 11; i++)
		EXPECT_EQ(Time_Display(pow(10.0, i)), results[i + 3]);
}

TEST(TestUtilities, TestPrintProgressBar)
{
	Print_Progress_Bar(0.4, 12);
}

TEST(TestUtilities, TestPrintBox)
{
	Print_Box("Hello", 1, 0);
}

// 2. Import and export data from files
TEST(TestUtilities, TestFileExists)
{
	// ARRANGE
	std::string file_name = "file_exists_test.txt";
	// ACT & ASSERT
	EXPECT_FALSE(File_Exists(file_name));
	std::ofstream f;
	f.open(file_name);
	f.close();
	EXPECT_TRUE(File_Exists(file_name));
	std::remove(file_name.c_str());
	EXPECT_FALSE(File_Exists(file_name));
}

TEST(TestUtilities, TestExportImportList)
{
	// ARRANGE
	double unit				 = 2.3;
	std::vector<double> list = {1.0 * unit, 2.0 * unit, 3.0 * unit, 4.0 * unit, 5.0 * unit};
	std::string file_name	 = "test_list.txt";
	// ACT
	Export_List(file_name, list, unit);
	std::vector<double> imported_list = Import_List(file_name, unit);
	// ASSERT
	ASSERT_EQ(imported_list.size(), list.size());
	for(unsigned int i = 0; i < list.size(); i++)
		ASSERT_DOUBLE_EQ(imported_list[i], list[i]);
	std::remove(file_name.c_str());
}

TEST(TestUtilities, TestExportImportTable)
{
	// ARRANGE
	std::vector<double> units			   = {2.3, 3.1};
	std::vector<std::vector<double>> table = {{11.0 * units[0], 12.0 * units[1]}, {21.0 * units[0], 22.0 * units[1]}};
	std::string file_name				   = "test_table.txt";
	// ACT
	Export_Table(file_name, table, units);
	std::vector<std::vector<double>> imported_table = Import_Table(file_name, units);
	// ASSERT
	ASSERT_EQ(imported_table.size(), table.size());
	for(unsigned int i = 0; i < table.size(); i++)
	{
		ASSERT_EQ(imported_table[i].size(), table[i].size());
		for(unsigned int j = 0; j < table[i].size(); j++)
			ASSERT_DOUBLE_EQ(imported_table[i][j], table[i][j]);
	}
	std::remove(file_name.c_str());
}

double func(double x)
{
	return x * x;
}
TEST(TestUtilities, TestExportImportFunction)
{
	// ARRANGE
	double xMin										= 1.0;
	double xMax										= 3.0;
	unsigned int steps								= 3;
	std::vector<std::vector<double>> correct_values = {{1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0}};
	std::string file_name							= "test_function.txt";
	// ACT
	Export_Function(file_name, func, xMin, xMax, steps);
	std::vector<std::vector<double>> imported_values = Import_Table(file_name);
	// ASSERT
	ASSERT_EQ(imported_values.size(), correct_values.size());
	for(unsigned int i = 0; i < correct_values.size(); i++)
	{
		ASSERT_EQ(imported_values[i].size(), correct_values[i].size());
		for(unsigned int j = 0; j < correct_values[i].size(); j++)
			ASSERT_DOUBLE_EQ(imported_values[i][j], correct_values[i][j]);
	}
	std::remove(file_name.c_str());
}

// 3. Create list with equi-distant numbers in log-space
TEST(TestUtilities, TestRangeMax)
{
	// ARRANGE
	unsigned int max = 13;
	// ACT & ASSERT
	ASSERT_EQ(Range(max).size(), max);
	for(unsigned int i = 0; i < max; i++)
		ASSERT_EQ(Range(max)[i], i);
}

TEST(TestUtilities, TestRangeMinMax)
{
	// ARRANGE
	unsigned int min = 6;
	unsigned int max = 13;
	// ACT & ASSERT
	ASSERT_EQ(Range(min, max).size(), max - min);
	for(unsigned int i = 0; i < Range(min, max).size(); i++)
		ASSERT_EQ(Range(min, max)[i], min + i);
}

TEST(TestUtilities, TestRangeMaxMin)
{
	// ARRANGE
	unsigned int min = 13;
	unsigned int max = 6;
	// ACT & ASSERT
	ASSERT_EQ(Range(min, max).size(), min - max);
	for(unsigned int i = 0; i < Range(min, max).size(); i++)
		ASSERT_EQ(Range(min, max)[i], min - i);
}

TEST(TestUtilities, TestLinearSpace)
{
	// ARRANGE
	double min						 = 1.0;
	double max						 = 10.0;
	unsigned int steps				 = 10;
	std::vector<double> list_correct = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	// ACT
	std::vector<double> list = Linear_Space(min, max, steps);
	// ASSERT
	ASSERT_EQ(list.size(), list_correct.size());
	for(unsigned int i = 0; i < list_correct.size(); i++)
		ASSERT_DOUBLE_EQ(list[i], list_correct[i]);
}

TEST(TestUtilities, TestLogSpace)
{
	// ARRANGE
	double min						 = 1.0e1;
	double max						 = 1.0e10;
	unsigned int steps				 = 10;
	std::vector<double> list_correct = {1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10};
	// ACT
	std::vector<double> list = Log_Space(min, max, steps);
	// ASSERT
	ASSERT_EQ(list.size(), list_correct.size());
	for(unsigned int i = 0; i < list_correct.size(); i++)
		ASSERT_NEAR(list[i], list_correct[i], 1.0e-3);
}

// 4. Dual stream class to write onto terminal and a log file simultaneously.
TEST(TestUtilities, TestLogger)
{
	// ARRANGE
	std::string file = "log.txt";
	Logger logger(file);
	std::vector<int> numbers = {17, 1, 9};
	// ACT
	for(auto& number : numbers)
		logger << number << std::endl;
	std::ifstream f(file);
	int number;
	std::vector<int> output;
	while(f >> number)
		output.push_back(number);
	f.close();
	// ASSERT
	EXPECT_TRUE(File_Exists(file));
	EXPECT_EQ(numbers.size(), output.size());
	for(unsigned int i = 0; i < numbers.size(); i++)
		EXPECT_EQ(numbers[i], output[i]);
}

// 5. Configuration class

// 6. Other utilities
TEST(TestUtilities, TestLocateClosestLocation)
{
	// ARRANGE
	std::vector<double> sorted_list	  = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	std::vector<double> unsorted_list = {3.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
	std::vector<double> targets		  = {-10.0, 0.4, 0.6, 0.9, 1.1, 1.55, 3.7, 4.0, 4.6, 8.1, 8.6, 9.0, 9.1, 10.9, 1e4};
	std::vector<unsigned int> result  = {0, 0, 1, 1, 1, 2, 4, 4, 5, 8, 9, 9, 9, 9, 9};
	// ACT & ASSERT
	for(int i = 0; i < targets.size(); i++)
		EXPECT_EQ(Locate_Closest_Location(sorted_list, targets[i]), result[i]);
	EXPECT_DEATH(Locate_Closest_Location(unsorted_list, 3.1), "");
}

TEST(TestUtilities, TestCheckForError)
{
	// ARRANGE
	std::string error_message = "Error message";
	double a				  = 1.0;
	double b				  = 2.0;
	// ACT & ASSERT
	EXPECT_DEATH(Check_For_Error(a == a, "function", error_message), error_message);
	EXPECT_NO_FATAL_FAILURE(Check_For_Error(a == b, "function", error_message));
}

TEST(TestUtilities, TestCheckForWarning)
{
	// ARRANGE
	std::string warning_message = "Warning message";
	double a					= 1.0;
	double b					= 2.0;
	// ACT & ASSERT
	EXPECT_NO_FATAL_FAILURE(Check_For_Warning(a == b, "function", warning_message));
	EXPECT_NO_FATAL_FAILURE(Check_For_Warning(a == a, "function", warning_message));
}