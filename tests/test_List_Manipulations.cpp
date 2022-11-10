#include "gtest/gtest.h"

#include "libphysica/List_Manipulations.hpp"

using namespace libphysica;

TEST(TestUtilities, TestListsEqual1)
{
	// ARRANGE
	std::vector<int> list1 = {1, 2, 3};
	std::vector<int> list2 = {4, 5, 6};
	std::vector<int> list3 = {1, 2, 3, 4};
	std::vector<int> list4 = {4, 5, 7};
	std::vector<int> list5 = {4, 5, 6};
	// ACT & ASSERT
	EXPECT_TRUE(Lists_Equal(list1, list1));
	EXPECT_FALSE(Lists_Equal(list1, list2));
	EXPECT_FALSE(Lists_Equal(list1, list3));
	EXPECT_FALSE(Lists_Equal(list1, list4));
	EXPECT_TRUE(Lists_Equal(list2, list5));
}

TEST(TestUtilities, TestListsEqual2)
{
	// ARRANGE
	std::vector<std::vector<double>> list1 = {{1, 2, 3}, {4, 5, 6}};
	std::vector<std::vector<double>> list2 = {{1, 2, 3, 4}, {4, 5, 6}};
	std::vector<std::vector<double>> list3 = {{1, 2, 13}, {4, 5, 6}};
	std::vector<std::vector<double>> list4 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	std::vector<std::vector<double>> list5 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	// ACT & ASSERT
	EXPECT_TRUE(Lists_Equal(list1, list1));
	EXPECT_FALSE(Lists_Equal(list1, list2));
	EXPECT_FALSE(Lists_Equal(list1, list3));
	EXPECT_FALSE(Lists_Equal(list1, list4));
	EXPECT_FALSE(Lists_Equal(list1, list5));
	EXPECT_TRUE(Lists_Equal(list4, list5));
}

TEST(TestUtilities, TestCombineLists)
{
	// ARRANGE
	std::vector<int> list1		   = {1, 2, 3};
	std::vector<int> list2		   = {4, 5, 6};
	std::vector<int> combined_list = {1, 2, 3, 4, 5, 6};
	// ACT & ASSERT
	ASSERT_TRUE(Lists_Equal(Combine_Lists(list1, list2), combined_list));
}

TEST(TestUtilities, TestTransposeLists)
{
	// ARRANGE
	std::vector<int> list1						  = {1, 2, 3};
	std::vector<int> list2						  = {4, 5, 6};
	std::vector<std::vector<int>> transposed_list = {{1, 4}, {2, 5}, {3, 6}};
	// ACT & ASSERT
	ASSERT_TRUE(Lists_Equal(Transpose_Lists(list1, list2), transposed_list));
}

TEST(TestUtilities, TestTransposeLists2)
{
	// ARRANGE
	std::vector<std::vector<int>> lists			  = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}};
	std::vector<std::vector<int>> transposed_list = {{1, 5, 9}, {2, 6, 10}, {3, 7, 11}, {4, 8, 12}};
	// ACT & ASSERT
	ASSERT_TRUE(Lists_Equal(Transpose_Lists(lists), transposed_list));
}

TEST(TestUtilities, TestSubList)
{
	// ARRANGE
	std::vector<int> list	 = {1, 2, 3, 4, 5, 6};
	std::vector<int> sublist = {3, 4, 5};
	// ACT & ASSERT
	ASSERT_TRUE(Lists_Equal(Sub_List(list, 2, 4), sublist));
}

TEST(TestUtilities, TestFlattenList)
{
	// ARRANGE
	std::vector<std::vector<std::string>> list = {{"a", "b", "c"}, {"d", "e", "f"}};
	std::vector<std::string> flattened_list	   = {"a", "b", "c", "d", "e", "f"};
	// ACT & ASSERT
	ASSERT_TRUE(Lists_Equal(Flatten_List(list), flattened_list));
}

TEST(TestUtilities, TestListContains)
{
	// ARRANGE
	std::vector<int> list = {1, 2, 3, 4, 5, 6};
	// ACT & ASSERT
	EXPECT_FALSE(List_Contains(list, 0));
	EXPECT_TRUE(List_Contains(list, 1));
	EXPECT_TRUE(List_Contains(list, 2));
	EXPECT_TRUE(List_Contains(list, 3));
	EXPECT_TRUE(List_Contains(list, 4));
	EXPECT_TRUE(List_Contains(list, 5));
	EXPECT_TRUE(List_Contains(list, 6));
	EXPECT_FALSE(List_Contains(list, 7));
}

TEST(TestUtilities, TestFindIndices)
{
	// ARRANGE
	std::vector<int> list = {1, 2, 3, 4, 3, 4, 5, 6, 6};
	// ACT & ASSERT
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 1), std::vector<int> {0}));
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 2), std::vector<int> {1}));
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 3), std::vector<int> {2, 4}));
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 4), std::vector<int> {3, 5}));
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 5), std::vector<int> {6}));
	EXPECT_TRUE(Lists_Equal(Find_Indices(list, 6), std::vector<int> {7, 8}));
}