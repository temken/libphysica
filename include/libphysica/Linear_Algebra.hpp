#ifndef __Linear_Algebra_hpp_
#define __Linear_Algebra_hpp_

#include <iostream>
#include <vector>

namespace libphysica
{
// 1. Vector class
class Vector
{
  private:
	std::vector<double> components;
	unsigned int dimension;

  public:
	// Constructors
	Vector();
	explicit Vector(unsigned int dim);
	Vector(unsigned int dim, double entry);
	Vector(std::vector<double> entries);
	Vector(const Vector& rhs);

	// Size functions
	unsigned int Size() const;
	void Resize(unsigned int dim);
	void Assign(unsigned int dim, double entry);

	// Products
	double Dot(const Vector& rhs) const;
	Vector Cross(const Vector& rhs) const;

	// Norm functions
	double Norm() const;
	void Normalize();
	Vector Normalized() const;

	// Overloading brackets
	double& operator[](const unsigned int i);
	const double& operator[](const unsigned int i) const;

	// Overloading operators
	Vector operator+(Vector v) const;
	Vector operator-(Vector v) const;
	double operator*(Vector v) const;
	Vector operator*(double s) const;
	Vector operator/(double s) const;
	Vector operator=(Vector v);
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);
};

Vector operator*(double s, const Vector& v);
std::ostream& operator<<(std::ostream& output, const Vector& v);

bool operator==(const Vector& v1, const Vector& v2);

// 2. Coordinate systems
extern Vector Spherical_Coordinates(double r, double theta, double phi);

// 3. Matrices
class Matrix
{
  private:
	std::vector<std::vector<double>> components;
	unsigned int rows, columns;

  public:
	// Constructors
	Matrix();
	Matrix(unsigned int dim_rows, unsigned int dim_columns);
	Matrix(unsigned int dim_rows, unsigned int dim_columns, double entry);
	Matrix(std::vector<std::vector<double>> entries);
	Matrix(const Matrix& rhs);
	Matrix(std::vector<double> diagonal_entries);
	Matrix(std::vector<std::vector<Matrix>> block_matrices);
	// Functions
	// Size and components
	unsigned int Rows() const;
	unsigned int Columns() const;
	void Resize(int row, int col);
	void Assign(int row, int col, double entry);
	void Delete_Row(unsigned int row);
	void Delete_Column(unsigned int column);
	Vector Return_Row(unsigned int row) const;
	Vector Return_Column(unsigned int column) const;

	//Binary operations
	Matrix Plus(const Matrix& M) const;
	Matrix Minus(const Matrix& M) const;
	Matrix Product(double s) const;
	Matrix Product(const Matrix& M) const;
	Vector Product(const Vector& v) const;
	Matrix Division(double s) const;

	//Matrix properties
	bool Square() const;
	bool Symmetric() const;
	bool Antisymmetric() const;
	bool Diagonal() const;
	bool Invertible() const;
	bool Orthogonal() const;

	//Matrix operations
	Matrix Transpose() const;
	Matrix Inverse() const;
	double Trace() const;
	double Determinant() const;

	Matrix Sub_Matrix(int row, int column) const;
	std::vector<double> Eigen_Values() const;	 //TODO
	std::vector<Vector> Eigen_Vectors() const;	 //TODO

	// Overloading brackets
	std::vector<double>& operator[](const unsigned int i);
	const std::vector<double>& operator[](const unsigned int i) const;
	// Overloading operators
	Matrix operator+(Matrix v);
	Matrix operator-(Matrix v);
	Matrix operator*(const Matrix& M);
	Vector operator*(const Vector& v_rhs);
	Matrix operator*(double s);
	Matrix operator/(double s);
	Matrix operator=(Matrix v);
	Matrix& operator+=(const Matrix& M);
	Matrix& operator-=(const Matrix& M);
};

// Overloading operators (Non-members)
Matrix operator*(double s, const Matrix& M);
Vector operator*(const Vector& v_left, const Matrix& M);
std::ostream& operator<<(std::ostream& output, const Matrix& v);
bool operator==(const Matrix& v1, const Matrix& v2);

extern Matrix Identity_Matrix(unsigned int dim);
extern Matrix Rotation_Matrix(double alpha, int dim, Vector axis = Vector({0, 0, 1}));
extern Matrix Outer_Vector_Product(const Vector& lhs, const Vector& rhs);

// Eigen systems
extern std::pair<Matrix, Matrix> QR_Decomposition(const Matrix& M);
extern std::vector<double> Eigenvalues(const Matrix& M);
extern std::vector<Vector> Eigenvectors(Matrix& M);
extern std::pair<std::vector<double>, std::vector<Vector>> Eigensystem(Matrix& M);

}	// namespace libphysica

#endif