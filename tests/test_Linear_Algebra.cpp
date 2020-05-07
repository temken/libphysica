#include "gtest/gtest.h"

#include <cmath>

#include "Linear_Algebra.hpp"

// 1. Vector class
	TEST(TestVector, TestSize)
	{
		// ARRANGE
		Vector v(4);
		unsigned int correct_dim = 4;
		// ACT & ASSERT
		ASSERT_EQ( v.Size() ,  correct_dim );
	}

	TEST(TestVector, TestResize)
	{
		// ARRANGE
		Vector v({1,0,0,0});
		unsigned int correct_result = 5;
		// ACT
		v.Resize(5);
		// ASSERT
		ASSERT_EQ( v.Size(), correct_result );
	}

	TEST(TestVector, TestAssign)
	{
		// ARRANGE
		Vector v({1,0,0,0});
		unsigned int dim = 10;
		double new_entry = 0.6;
		// ACT
		v.Assign(dim, new_entry);
		// ASSERT
		ASSERT_EQ( v.Size(), dim );
		for(unsigned int i = 0; i < dim; i++)
			ASSERT_DOUBLE_EQ( v[i], new_entry );
	}

	TEST(TestVector, TestNorm)
	{
		// ARRANGE
		Vector v(3 , 2.0);
		double correct_norm = sqrt(12.0);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( v.Norm(), correct_norm );
	}

	TEST(TestVector, TestDot)
	{
		// ARRANGE
		Vector v1({1,2,3});
		Vector v2(3, 2.0);
		double correct_result = 12.0;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( v1.Dot(v2) ,  correct_result );
	}

	TEST(TestVector, TestCross)
	{
		// ARRANGE
		Vector v1(3 , 2.0);
		Vector v2(v1);
		Vector correct_result(3, 0.0);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( v1.Cross(v2)[0], correct_result[0] );
		ASSERT_DOUBLE_EQ( v1.Cross(v2)[1], correct_result[1] );
		ASSERT_DOUBLE_EQ( v1.Cross(v2)[2], correct_result[2] );
	}


	TEST(TestVector, TestNormalize)
	{
		// ARRANGE
		Vector v(3 , 2.0);
		double correct_norm = 1.0;
		// ACT
		v.Normalize();
		// ASSERT
		ASSERT_DOUBLE_EQ( v.Norm(), correct_norm );
	}

	TEST(TestVector, TestNormalized)
	{
		// ARRANGE
		Vector v({100,100,100});
		double correct_norm = 1.0;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( v.Normalized().Norm(), correct_norm );
	}

	TEST(TestVector, TestBrackets)
	{
		// ARRANGE
		Vector v1({2.5,5.0,7.5});
		const Vector v2({2.5,5.0,7.5});
		double correct_entry_0 = 2.5;
		double correct_entry_1 = 5.0;
		double correct_entry_2 = 7.5;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( v1[0], correct_entry_0 );
		ASSERT_DOUBLE_EQ( v1[1], correct_entry_1 );
		ASSERT_DOUBLE_EQ( v1[2], correct_entry_2 );
		
		ASSERT_DOUBLE_EQ( v2[0], correct_entry_0 );
		ASSERT_DOUBLE_EQ( v2[1], correct_entry_1 );
		ASSERT_DOUBLE_EQ( v2[2], correct_entry_2 );
	}

	TEST(TestVector, TestOperators)
	{
		// ARRANGE
		Vector v1({2.5,5.0,7.5});
		Vector v2(3,1.0);
		Vector v3 = v1;
		double s = 2.0;

		Vector v1_times_s({5.0,10.0,15.0});
		Vector v2_times_s({2.0,2.0,2.0});
		Vector v1_by_s({1.25,2.5,3.75});
		Vector v1_plus_v2({3.5,6.0,8.5});
		Vector v1_minus_v2({1.5,4.0,6.5});
		double v1_times_v2 = 15.0;

		// ACT & ASSERT
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( (v1 + v2)[i], v1_plus_v2[i]);
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( (v1 - v2)[i], v1_minus_v2[i]);
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( (s * v1)[i], v1_times_s[i]);
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( (v2 * s)[i], v2_times_s[i]);
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( (v1 / s)[i], v1_by_s[i]);
		ASSERT_DOUBLE_EQ( v1 * v2, v1_times_v2);
		v3 += v2;
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( v3[i], v1_plus_v2[i]);
		v3 -= v1;
		for(unsigned int i = 0; i < 3; i++)
			ASSERT_DOUBLE_EQ( v3[i], v2[i]);
		ASSERT_TRUE( v2 == v3 );
	}

// 2. Coordinate systems
	TEST(TestCoordinateSystems, TestSphericalCoordinates)
	{
		// ARRANGE
		double r = 3.14;
		double theta = 0.4;
		double phi = 0.5;
		Vector R = Spherical_Coordinates(r, theta, phi);
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ( R.Norm(), r );
		ASSERT_DOUBLE_EQ( acos(R[2] / R.Norm()), theta);
		ASSERT_DOUBLE_EQ( atan2(R[1], R[0]), phi);
	}

// 3. Matrices
	TEST(TestMatrix, TestRowsColumns)
	{
		// ARRANGE
		Matrix M(2,3);
		unsigned int rows = 2;
		unsigned int columns = 3;
		// ACT & ASSERT
		ASSERT_EQ(M.Rows(), rows);
		ASSERT_EQ(M.Columns(), columns);
	}

	TEST(TestMatrix, TestResize)
	{
		// ARRANGE
		Matrix M(2,3,0.0);
		M.Resize(5,6);
		unsigned int rows = 5;
		unsigned int columns = 6;
		// ACT & ASSERT
		ASSERT_EQ(M.Rows(), rows);
		ASSERT_EQ(M.Columns(), columns);
	}

	TEST(TestMatrix, TestAssign)
	{
		// ARRANGE
		Matrix M({{1,0},{0,1}});
		unsigned int rows = 5;
		unsigned int columns = 6;
		double entry = 2.3;
		M.Assign(rows, columns, entry);
		// ACT & ASSERT
		ASSERT_EQ(M.Rows(), rows);
		ASSERT_EQ(M.Columns(), columns);
		for(unsigned int i = 0; i < rows; i++)
			for(unsigned int j = 0; j < columns; j++)
				ASSERT_DOUBLE_EQ( M[i][j], entry);
	}

	TEST(TestMatrix, TestDeleteRow)
	{
		// ARRANGE
		Matrix M({{1,1,1},{2,2,2},{3,3,3}});
		Matrix correct_result({{1,1,1},{3,3,3}});
		unsigned int rows = 2;
		unsigned int columns = 3;
		// ACT
		M.Delete_Row(1);
		// ASSERT
		ASSERT_EQ(M.Rows(), rows);
		ASSERT_EQ(M.Columns(), columns);
		for(unsigned int i = 0; i < rows; i++)
			for(unsigned int j = 0; j < columns; j++)
				ASSERT_DOUBLE_EQ( M[i][j], correct_result[i][j]);
	}

	TEST(TestMatrix, TestDeleteColumn)
	{
		// ARRANGE
		Matrix M({{1,1,1},{2,2,2},{3,3,3}});
		Matrix correct_result({{1,1},{2,2},{3,3}});
		unsigned int rows = 3;
		unsigned int columns = 2;
		// ACT
		M.Delete_Column(1);
		// ASSERT
		ASSERT_EQ(M.Rows(), rows);
		ASSERT_EQ(M.Columns(), columns);
		for(unsigned int i = 0; i < rows; i++)
			for(unsigned int j = 0; j < columns; j++)
				ASSERT_DOUBLE_EQ( M[i][j], correct_result[i][j]);
	}

	TEST(TestMatrix, TestPlus)
	{
		// ARRANGE
		Matrix M1({{1,1,1},{2,2,2},{3,3,3}});
		Matrix M2 = Unit_Matrix(3);
		Matrix M1_plus_M2({{2,1,1},{2,3,2},{3,3,4}});
		unsigned int rows = 3;
		unsigned int columns = 3;
		// ACT
		Matrix M3 = M1.Plus(M2);
		// ASSERT
		ASSERT_EQ(M3.Rows(), rows);
		ASSERT_EQ(M3.Columns(), columns);
		for(unsigned int i = 0; i < M3.Rows(); i++)
			for(unsigned int j = 0; j < M3.Columns(); j++)
				ASSERT_DOUBLE_EQ(M3[i][j], M1_plus_M2[i][j]);
	}

	TEST(TestMatrix, TestMinus)
	{
		// ARRANGE
		Matrix M1({{1,1,1},{2,2,2},{3,3,3}});
		Matrix M2 = Unit_Matrix(3);
		Matrix M1_minus_M2({{0,1,1},{2,1,2},{3,3,2}});
		unsigned int rows = 3;
		unsigned int columns = 3;
		// ACT
		Matrix M3 = M1.Minus(M2);
		// ASSERT
		ASSERT_EQ(M3.Rows(), rows);
		ASSERT_EQ(M3.Columns(), columns);
		for(unsigned int i = 0; i < M3.Rows(); i++)
			for(unsigned int j = 0; j < M3.Columns(); j++)
				ASSERT_DOUBLE_EQ(M3[i][j], M1_minus_M2[i][j]);
	}

	TEST(TestMatrix, TestProductScalar)
	{
		// ARRANGE
		Matrix M({{1,1,1},{2,2,2},{3,3,3}});
		double s = 2.5;
		Matrix correct_result({{2.5,2.5,2.5},{5.0,5.0,5.0},{7.5,7.5,7.5}});
		// ACT
		Matrix M2 = M.Product(s);
		// ASSERT
		for(unsigned int i = 0; i < M2.Rows(); i++)
			for(unsigned int j = 0; j < M2.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], correct_result[i][j]);
	}


	TEST(TestMatrix, TestProductMatrix)
	{
		// ARRANGE
		Matrix M1({{1,1,1},{2,2,2},{3,3,3}});
		Matrix M2({{4,4,4},{5,5,5},{6,6,6}});
		Matrix M1_times_M2({{15, 15, 15}, {30, 30, 30}, {45, 45, 45}});
		// ACT
		Matrix M3 = M1.Product(M2);
		// ASSERT
		for(unsigned int i = 0; i < M1.Rows(); i++)
			for(unsigned int j = 0; j < M1.Columns(); j++)
				ASSERT_DOUBLE_EQ(M3[i][j], M1_times_M2[i][j]);
	}

	TEST(TestMatrix, TestProductVector)
	{
		// ARRANGE
		Matrix M({{1,1,1},{2,2,2},{3,3,3}});
		Vector v1({4,5,6});
		Vector M_times_v1({15, 30, 45});
		// ACT
		Vector v2 = M.Product(v1);
		// ASSERT
		for(unsigned int i = 0; i < v1.Size(); i++)
			ASSERT_DOUBLE_EQ(v2[i], M_times_v1[i]);
	}

	TEST(TestMatrix, TestDivision)
	{
		// ARRANGE
		Matrix M1({{1,1,1},{2,2,2},{3,3,3}});
		double s = 2.0;
		Matrix M1_by_s({{0.5,0.5,0.5},{1.0,1.0,1.0},{1.5,1.5,1.5}});
		// ACT
		Matrix M2 = M1.Division(s);
		// ASSERT
		for(unsigned int i = 0; i < M1.Rows(); i++)
			for(unsigned int j = 0; j < M1.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], M1_by_s[i][j]);
	}

	TEST(TestMatrix, TestSquare)
	{
		// ARRANGE
		Matrix M1 = Unit_Matrix(3);
		Matrix M2(3,2);
		// ACT & ASSERT
		ASSERT_TRUE(M1.Square());
		ASSERT_FALSE(M2.Square());
	}

	TEST(TestMatrix, TestSymmetric)
	{
		// ARRANGE
		Matrix M1({{1,2,3},{2,4,5},{3,5,6}});
		Matrix M2({{1,2,3},{4,5,6},{7,8,9}});
		// ACT & ASSERT
		ASSERT_TRUE(M1.Symmetric());
		ASSERT_FALSE(M2.Symmetric());
	}

	TEST(TestMatrix, TestAntisymmetric)
	{
		// ARRANGE
		Matrix M1({{0,2,3},{-2,0,5},{-3,-5,0}});
		Matrix M2({{1,2,3},{2,4,5},{3,5,6}});
		// ACT & ASSERT
		ASSERT_TRUE(M1.Antisymmetric());
		ASSERT_FALSE(M2.Antisymmetric());
	}

	TEST(TestMatrix, TestDiagonal)
	{
		// ARRANGE
		Matrix M1({{3,0,0},{0,2,0},{0,0,-3.3}});
		Matrix M2({{1,2,3},{2,4,5},{3,5,6}});
		// ACT & ASSERT
		ASSERT_TRUE(M1.Diagonal());
		ASSERT_FALSE(M2.Diagonal());
	}

	TEST(TestMatrix, TestInvertible)
	{
		// ARRANGE
		Matrix M1({{3,2,1},{4,2,1},{3,10,-3.3}});
		Matrix M2({{3,0,0},{0,2,0},{0,0,0}});
		// ACT & ASSERT
		ASSERT_TRUE(M1.Invertible());
		ASSERT_FALSE(M2.Invertible());
	}

	// TEST(TestMatrix, TestOrthogonal)
	// {
	// 	// ARRANGE
	// 	Matrix M1({{1.0/sqrt(5.0),2.0/sqrt(5.0)},{2.0/sqrt(5.0),-1.0/sqrt(5.0)}});
	// 	Matrix M2({{3,2},{0,2}});
	// 	// ACT & ASSERT
	// 	ASSERT_TRUE(M1.Orthogonal());
	// 	ASSERT_FALSE(M2.Orthogonal());
	// }

	TEST(TestMatrix, TestTranspose)
	{
		// ARRANGE
		Matrix M1({{1,2,3},{4,5,6},{7,8,9}});
		Matrix M_Transposed({{1,4,7},{2,5,8},{3,6,9}});
		// ACT 
		Matrix M2 = M1.Transpose();
		// ASSERT
		for(unsigned int i = 0; i < M1.Rows(); i++)
			for(unsigned int j = 0; j < M1.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], M_Transposed[i][j]);
	}

	TEST(TestMatrix, TestInverse)
	{
		// ARRANGE
		Matrix M1({{1,0,0},{0,2,0},{0,0,3}});
		Matrix M_Inverse({{1,0,0},{0,0.5,0},{0,0,1.0/3.0}});
		// ACT 
		Matrix M2 = M1.Inverse();
		// ASSERT
		for(unsigned int i = 0; i < M1.Rows(); i++)
			for(unsigned int j = 0; j < M1.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], M_Inverse[i][j]);
	}

	TEST(TestMatrix, TestTrace)
	{
		// ARRANGE
		Matrix M1(4,4,2.0);
		double correct_trace = 8.0;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ(M1.Trace(), correct_trace);
	}

	TEST(TestMatrix, TestDeterminant)
	{
		// ARRANGE
		Matrix M1({{3,2,1},{4,2,1},{3,10,-3.3}});
		double correct_determinant = 16.6;
		// ACT & ASSERT
		ASSERT_DOUBLE_EQ(M1.Determinant(), correct_determinant);
	}

	TEST(TestMatrix, TestSubMatrix)
	{
		// ARRANGE
		Matrix M1({{11,12,13},{21,22,23},{31,32,33}});
		Matrix correct_submatrix({{11,13},{31,33}});
		// ACT
		Matrix M2 = M1.Sub_Matrix(1,1);
		// ASSERT
		for(unsigned int i = 0; i < M2.Rows(); i++)
			for(unsigned int j = 0; j < M2.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], correct_submatrix[i][j]);
	}

	TEST(TestMatrix, TestBrackets)
	{
		// ARRANGE
		std::vector<std::vector<double>> entries = {{11,12,13},{21,22,23},{31,32,33}};
		Matrix M(entries);
		// ACT & ASSERT
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M[i][j], entries[i][j]);
	}

	TEST(TestMatrix, TestOperators)
	{
		// ARRANGE
		Matrix M({{11,12,13},{21,22,23},{31,32,33}});
		Vector v({1,2,3});
		double s = 1.5;

		Matrix M_plus_M({{22,24,26},{42,44,46},{62,64,66}});
		Matrix M_minus_M(3,3,0.0);
		Matrix M_times_M({{776, 812, 848}, {1406, 1472, 1538}, {2036, 2132, 2228}});
		Vector M_times_v({74, 134, 194});
		Matrix M_times_s({{16.5, 18., 19.5}, {31.5, 33., 34.5}, {46.5, 48., 49.5}});
		Matrix M_by_s({{22.0/3.0, 8.0, 26.0/3.0}, {14.0, 44.0/3.0, 46.0/3.0}, {62.0/3.0, 64.0/3.0, 22.0}});
		// ACT
		Matrix M2 = M + M;
		Matrix M3 = M - M;
		Matrix M4 = M * M;
		Vector v1 = M * v;
		Matrix M5 = s * M;
		Matrix M6 = M * s;
		Matrix M7 = M / s;

		Matrix M_copy = M;
		// ASSERT
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], M_plus_M[i][j]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M3[i][j], M_minus_M[i][j]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M4[i][j], M_times_M[i][j]);
		for(unsigned int i = 0; i < v.Size(); i++)
			ASSERT_DOUBLE_EQ(v1[i], M_times_v[i]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M5[i][j], M_times_s[i][j]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M6[i][j], M_times_s[i][j]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M7[i][j], M_by_s[i][j]);
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M_copy[i][j], M[i][j]);
		M_copy += M;
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M_copy[i][j], M_plus_M[i][j]);
		M_copy -= M+M;
		for(unsigned int i = 0; i < M.Rows(); i++)
			for(unsigned int j = 0; j < M.Columns(); j++)
				ASSERT_DOUBLE_EQ(M_copy[i][j], M_minus_M[i][j]);
	}

	TEST(TestMatrix, TestUnitMatrix)
	{
		// ARRANGE

		Matrix unit_matrix_3({{1,0,0},{0,1,0},{0,0,1}});
		Matrix unit_matrix_4({{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}});
		// ACT
		Matrix M1;
		Matrix M2 = Unit_Matrix(4);
		// ASSERT
		for(unsigned int i = 0; i < unit_matrix_3.Rows(); i++)
			for(unsigned int j = 0; j < unit_matrix_3.Columns(); j++)
				ASSERT_DOUBLE_EQ(M1[i][j], unit_matrix_3[i][j]);
		for(unsigned int i = 0; i < unit_matrix_4.Rows(); i++)
			for(unsigned int j = 0; j < unit_matrix_4.Columns(); j++)
				ASSERT_DOUBLE_EQ(M2[i][j], unit_matrix_4[i][j]);
	}

	TEST(TestMatrix, TestRotationMatrix)
	{
		// ARRANGE
		Matrix rotation_matrix({{1.0/sqrt(2.0),-1.0/sqrt(2.0)},{1.0/sqrt(2.0),1.0/sqrt(2.0)}});
		// ACT
		Matrix M = Rotation_Matrix(M_PI/4.0, 2);
		// ASSERT
		for(unsigned int i = 0; i < rotation_matrix.Rows(); i++)
			for(unsigned int j = 0; j < rotation_matrix.Columns(); j++)
				ASSERT_DOUBLE_EQ(M[i][j], rotation_matrix[i][j]);
	}