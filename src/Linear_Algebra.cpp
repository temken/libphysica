#include "Linear_Algebra.hpp"

#include <cmath>

#include "Numerics.hpp"

namespace libphysica
{

// 1. Vector functions.
// Constructors
Vector::Vector()
: Vector(3)
{
}

Vector::Vector(unsigned int dim)
: Vector(dim, 0.0)
{
}

Vector::Vector(unsigned int dim, double entry)
: dimension(dim)
{
	components = std::vector<double>(dim, entry);
}

Vector::Vector(std::vector<double> entries)
: components(entries)
{
	dimension = entries.size();
}

Vector::Vector(const Vector& rhs)
: components(rhs.components), dimension(rhs.dimension)
{
}

// Functions
unsigned int Vector::Size() const
{
	return dimension;
}

void Vector::Resize(unsigned int dim)
{
	dimension = dim;
	components.resize(dim);
}

void Vector::Assign(unsigned int dim, double entry)
{
	dimension = dim;
	components.assign(dim, entry);
}

double Vector::Norm() const
{
	return sqrt(Dot(*this));
}

double Vector::Dot(const Vector& rhs) const
{
	if(dimension != rhs.Size())
	{
		std::cerr << "Error in Vector::Dot(): Vectors are of differing dimensions (" << dimension << " and " << rhs.Size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		double result = 0;
		for(unsigned int i = 0; i < dimension; i++)
		{
			result += components[i] * rhs[i];
		}
		return result;
	}
}

Vector Vector::Cross(const Vector& rhs) const
{
	if(dimension != 3)
	{
		std::cerr << "Error in Vector Vector::Cross(): Cross product only defined for 3 dimensions, not " << dimension << "." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> result(3);
		result[0] = components[1] * rhs[2] - components[2] * rhs[1];
		result[1] = components[2] * rhs[0] - components[0] * rhs[2];
		result[2] = components[0] * rhs[1] - components[1] * rhs[0];
		return Vector(result);
	}
}

void Vector::Normalize()
{
	double norm = Norm();
	for(unsigned int i = 0; i < dimension; i++)
		components[i] = components[i] / norm;
}

Vector Vector::Normalized() const
{
	double norm = Norm();
	std::vector<double> new_components(dimension);
	for(unsigned int i = 0; i < dimension; i++)
		new_components[i] = (double) components[i] / norm;
	return Vector(new_components);
}

// Brackets
double& Vector::operator[](const unsigned int i)
{
	if(i < 0 || i >= dimension)
	{
		std::cerr << "Error in Vector::operator[]: Index " << i << " out of bound [" << 0 << "," << dimension - 1 << "]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return components[i];
}
const double& Vector::operator[](const unsigned int i) const
{
	if(i < 0 || i >= dimension)
	{
		std::cerr << "Error in Vector::operator[]: Index " << i << " out of bound [" << 0 << "," << dimension - 1 << "]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return components[i];
}

// Member operators
Vector Vector::operator+(Vector v)
{
	if(dimension != v.dimension)
	{
		std::cerr << "Error in Vector::operator+(Vector v): Two vectors of dimensions " << dimension << " and " << v.dimension << " can not be added." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> new_components(dimension);
		for(unsigned int i = 0; i < dimension; i++)
			new_components[i] = components[i] + v[i];
		return Vector(new_components);
	}
}

Vector Vector::operator-(Vector v)
{
	if(dimension != v.dimension)
	{
		std::cerr << "Error in Vector::operator-(Vector v): Two vectors of dimensions " << dimension << " and " << v.dimension << " can not be subtracted." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> new_components(dimension);
		for(unsigned int i = 0; i < dimension; i++)
			new_components[i] = components[i] - v[i];
		return Vector(new_components);
	}
}

double Vector::operator*(Vector v)
{
	return Dot(v);
}

Vector Vector::operator*(double s)
{
	std::vector<double> new_components(dimension);
	for(unsigned int i = 0; i < dimension; i++)
		new_components[i] = components[i] * s;
	return Vector(new_components);
}

Vector Vector::operator/(double s)
{
	std::vector<double> new_components(dimension);
	for(unsigned int i = 0; i < dimension; i++)
		new_components[i] = components[i] / s;
	return Vector(new_components);
}

Vector Vector::operator=(Vector v)
{
	components = v.components;
	dimension  = v.dimension;
	return *this;
}

Vector& Vector::operator+=(const Vector& v)
{
	for(unsigned int i = 0; i < dimension; i++)
		components[i] += v[i];
	return *this;
}

Vector& Vector::operator-=(const Vector& v)
{
	for(unsigned int i = 0; i < dimension; i++)
		components[i] -= v[i];
	return *this;
}

// Non-member operators
Vector operator*(double s, const Vector& v)
{
	std::vector<double> new_components(v.Size());
	for(unsigned int i = 0; i < v.Size(); i++)
		new_components[i] = v[i] * s;
	return Vector(new_components);
}

std::ostream& operator<<(std::ostream& output, const Vector& v)
{
	output << "(";
	for(unsigned int i = 0; i < v.Size(); i++)
	{
		output << v[i];
		if(i < (v.Size() - 1))
			output << " , ";
	}
	output << ")";
	return output;
}

bool operator==(const Vector& v1, const Vector& v2)
{
	if(v1.Size() != v2.Size())
		return false;
	else
	{
		bool result = true;
		for(unsigned int i = 0; i < v1.Size(); i++)
		{
			result *= (v1[i] == v2[i]);
			if(!result)
				break;
		}
		return result;
	}
}

// 2. Coordinates
Vector Spherical_Coordinates(double r, double theta, double phi)
{
	std::vector<double> comp = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
	return Vector(comp);
}

// 3. Matrix functions.
// Constructors
Matrix::Matrix()
: rows(3), columns(3)
{
	components.push_back(std::vector<double>{1.0, 0, 0});
	components.push_back(std::vector<double>{0, 1.0, 0});
	components.push_back(std::vector<double>{0, 0, 1.0});
}

Matrix::Matrix(unsigned int dim_rows, unsigned int dim_columns)
: Matrix(dim_rows, dim_columns, 0.0)
{
}

Matrix::Matrix(unsigned int dim_rows, unsigned int dim_columns, double entry)
: rows(dim_rows), columns(dim_columns)
{
	for(unsigned int i = 0; i < rows; i++)
	{
		std::vector<double> aux(dim_columns, entry);
		components.push_back(aux);
	}
}

Matrix::Matrix(std::vector<std::vector<double>> entries)
: components(entries), rows(entries.size()), columns(entries[0].size())
{
	for(unsigned int i = 0; i < rows; i++)
	{
		if(entries[i].size() != columns)
		{
			std::cerr << "Error in Matrix::Matrix(): The input vector<vector<double>> has irregular shape." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	components = entries;
}

Matrix::Matrix(const Matrix& rhs)
: components(rhs.components), rows(rhs.rows), columns(rhs.columns)
{
}

Matrix::Matrix(std::vector<double> diagonal_entries)
{
	rows	   = diagonal_entries.size();
	columns	   = diagonal_entries.size();
	components = std::vector<std::vector<double>>(rows, std::vector<double>(columns, 0.0));
	for(unsigned int i = 0; i < rows; i++)
		components[i][i] = diagonal_entries[i];
}

// Functions
// Size and components
unsigned int Matrix::Rows() const
{
	return rows;
}

unsigned int Matrix::Columns() const
{
	return columns;
}

void Matrix::Resize(int row, int col)
{
	rows	= row;
	columns = col;
	components.resize(row);
	for(unsigned int i = 0; i < rows; i++)
	{
		components[i].resize(col);
	}
}

void Matrix::Assign(int row, int col, double entry)
{
	rows	= row;
	columns = col;
	components.resize(row);
	for(unsigned int i = 0; i < rows; i++)
	{
		components[i].assign(col, entry);
	}
}

void Matrix::Delete_Row(unsigned int row)
{
	if(row < 0 || row >= rows)
	{
		std::cerr << "Error in Matrix::Delete_Row(int): Row " << row << " of a " << rows << "x" << columns << " matrix does not exist.\n(Note: The index domain starts at 0)" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		rows--;
		components.erase(components.begin() + row);
	}
}

void Matrix::Delete_Column(unsigned int column)
{
	if(column < 0 || column >= columns)
	{
		std::cerr << "Error in Matrix::Delete_Column(int): Column " << column << " of a " << rows << "x" << columns << " matrix does not exist.\n(Note: The index domain starts at 0)" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{

		columns--;
		for(unsigned int i = 0; i < rows; i++)
			components[i].erase(components[i].begin() + column);
	}
}
//Binary operations
Matrix Matrix::Plus(const Matrix& M) const
{
	if(rows != M.Columns() || columns != M.Rows())
	{
		std::cerr << "Error in Matrix::Plus(const Matrix&): The dimensions of the two matrices to be added do not match: (" << rows << "x" << columns << ") + (" << M.Rows() << "x" << M.Columns() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<std::vector<double>> result_components(rows, std::vector<double>(columns, 0.0));
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = 0; j < columns; j++)
			{
				result_components[i][j] = components[i][j] + M[i][j];
			}
		}
		return Matrix(result_components);
	}
}

Matrix Matrix::Minus(const Matrix& M) const
{
	if(rows != M.Columns() || columns != M.Rows())
	{
		std::cerr << "Error in Matrix::Minus(const Matrix&): The dimensions of the two matrices to be subtracted do not match: (" << rows << "x" << columns << ") - (" << M.Rows() << "x" << M.Columns() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<std::vector<double>> result_components(rows, std::vector<double>(columns, 0.0));
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = 0; j < columns; j++)
			{
				result_components[i][j] = components[i][j] - M[i][j];
			}
		}
		return Matrix(result_components);
	}
}

Matrix Matrix::Product(double s) const
{
	std::vector<std::vector<double>> result_components(rows, std::vector<double>(columns, 0.0));
	for(unsigned int i = 0; i < rows; i++)
	{
		for(unsigned int j = 0; j < columns; j++)
		{
			result_components[i][j] = s * components[i][j];
		}
	}
	return Matrix(result_components);
}

Matrix Matrix::Product(const Matrix& M) const
{
	if(columns != M.Rows())
	{
		std::cerr << "Error in Matrix::Product(const Matrix &M): # of columns of the lhs does not match # of rows of the rhs matrix: (" << rows << "x" << columns << ")*(" << M.Rows() << "x" << M.Columns() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		Matrix result(rows, M.Columns(), 0.0);

		for(unsigned int i = 0; i < result.Rows(); i++)
		{
			for(unsigned int j = 0; j < result.Columns(); j++)
			{
				for(unsigned int k = 0; k < columns; k++)
					result[i][j] += components[i][k] * M[k][j];
			}
		}
		return result;
	}
	return Matrix();
}

Vector Matrix::Product(const Vector& v_rhs) const
{
	if(v_rhs.Size() != columns)
	{
		std::cerr << "Error in Matrix::Product(const Vector&): The dimensions of the matrix and the vector do not match: (" << rows << "x" << columns << ")*(" << v_rhs.Size() << ")." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> result_components(rows, 0.0);
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = 0; j < columns; j++)
				result_components[i] += components[i][j] * v_rhs[j];
		}
		return Vector(result_components);
	}
}

Matrix Matrix::Division(double s) const
{
	std::vector<std::vector<double>> result_components(rows, std::vector<double>(columns, 0.0));

	for(unsigned int i = 0; i < rows; i++)
	{
		for(unsigned int j = 0; j < columns; j++)
		{
			result_components[i][j] = components[i][j] / s;
		}
	}
	return Matrix(result_components);
}

//Matrix properties
bool Matrix::Square() const
{
	if(rows == columns)
		return true;
	else
		return false;
}

bool Matrix::Symmetric() const
{
	if(!Square())
		return false;
	else
	{
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = i; j < columns; j++)
			{
				if(components[i][j] != components[j][i])
					return false;
			}
		}
		return true;
	}
}

bool Matrix::Antisymmetric() const
{
	if(!Square())
		return false;
	else
	{
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = i; j < columns; j++)
			{
				if(components[i][j] != -1.0 * components[j][i])
					return false;
			}
		}
		return true;
	}
}

bool Matrix::Diagonal() const
{
	if(!Square())
		return false;
	else
	{
		for(unsigned int i = 0; i < rows; i++)
		{
			for(unsigned int j = 0; j < columns; j++)
			{
				if(i != j && components[i][j] != 0.0)
					return false;
			}
		}
		return true;
	}
}

bool Matrix::Invertible() const
{
	if(!Square())
		return false;
	else
	{
		if(Determinant() != 0.0)
			return true;
		else
			return false;
	}
}

bool Matrix::Orthogonal() const
{
	if(!Invertible())
		return false;
	else
	{
		Matrix MT	= Transpose();
		Matrix Minv = Inverse();
		if(MT == Minv)
			return true;
		else
			return false;
	}
}

//Matrix operations
Matrix Matrix::Transpose() const
{
	unsigned int rows_new = columns;
	unsigned int col_new  = rows;
	std::vector<std::vector<double>> result_components(rows_new, std::vector<double>(col_new, 0.0));
	for(unsigned int i = 0; i < rows; i++)
	{
		for(unsigned int j = 0; j < columns; j++)
			result_components[j][i] = components[i][j];
	}
	return Matrix(result_components);
}

Matrix Matrix::Inverse() const
{
	if(!Square())
	{
		std::cerr << "Error in Matrix::Inverse(): Inverse matrix only defined for square matrices." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(!Invertible())
	{
		std::cerr << "Error in Matrix::Inverse(): Matrix is not invertible." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	// Gauss Jordan inversion
	else
	{
		// Augment the matrix with an identity matrix
		unsigned int N = rows;
		Matrix A(N, 2.0 * N, 0.0);

		for(unsigned int i = 0; i < N; i++)
		{
			for(unsigned int j = 0; j < N; j++)
			{
				A[i][j] = components[i][j];
				if(i == j)
					A[i][j + N] = 1.0;
				else
					A[i][j + N] = 0.0;
			}
		}
		// Gauss Jordan elimination
		for(unsigned int i = 0; i < N; i++)
		{
			if(A[i][i] == 0)
			{
				std::cerr << "Error in Matrix::Inverse(): Diagonal element is zero." << std::endl;
				std::exit(EXIT_FAILURE);
			}
			for(unsigned int j = 0; j < N; j++)
			{
				if(i != j)
				{
					double ratio = A[j][i] / A[i][i];
					for(unsigned int k = 0; k < A.Columns(); k++)
					{
						A[j][k] = A[j][k] - ratio * A[i][k];
					}
				}
			}
		}

		// Change diagonal entries to 1
		for(unsigned int i = 0; i < N; i++)
		{
			for(unsigned int j = N; j < 2 * N; j++)
			{
				A[i][j] = A[i][j] / A[i][i];
			}
			// A[i][i] /= A[i][i]; //not necessary
		}
		// Remove the first N columns.
		for(unsigned int i = 0; i < N; i++)
			A.Delete_Column(0);
		return A;
	}
}

double Matrix::Trace() const
{
	if(rows != columns)
	{
		std::cerr << "Error in Matrix::Trace(): Trace only defined for square matrices." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		double tr = 0.0;
		for(unsigned int i = 0; i < rows; i++)
			tr += components[i][i];
		return tr;
	}
}

double Matrix::Determinant() const
{
	if(!Square())
	{
		std::cerr << "Error in Matrix::Determinant(): Determinant only defined for square matrices." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else if(rows == 1)
	{
		return components[0][0];
	}
	else if(rows == 2)
	{
		return components[0][0] * components[1][1] - components[0][1] * components[1][0];
	}
	else
	{
		std::vector<double> factors(rows, 0.0);
		std::vector<Matrix> sub_matrices = {};
		for(unsigned int j = 0; j < columns; j++)
		{
			double sign = (j % 2 == 0) ? +1.0 : -1.0;
			factors[j]	= sign * components[0][j];
			sub_matrices.push_back(Sub_Matrix(0, j));
		}
		double det = 0.0;
		for(unsigned int j = 0; j < columns; j++)
			det += factors[j] * sub_matrices[j].Determinant();
		return det;
	}
}

Matrix Matrix::Sub_Matrix(int row, int column) const
{
	Matrix M(components);
	M.Delete_Row(row);
	M.Delete_Column(column);
	return M;
}

std::vector<double> Matrix::Eigen_Values() const
//TO DO
{
	return std::vector<double>(1);
}

std::vector<Vector> Matrix::Eigen_Vectors() const
//TO DO
{
	return std::vector<Vector>(1);
}

// Overloading brackets
std::vector<double>& Matrix::operator[](const unsigned int i)
{
	if(i < 0 || i > (rows - 1))
	{
		std::cerr << "Error in Matrix::operator[](): Index i=" << i << " is out of bound [" << 0 << "," << (rows - 1) << "]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return components[i];
}
const std::vector<double>& Matrix::operator[](const unsigned int i) const
{
	if(i < 0 || i > (rows - 1))
	{
		std::cerr << "Error in Matrix::operator[](): Index i=" << i << " is out of bound [" << 0 << "," << (rows - 1) << "]." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
		return components[i];
}

// Overloading operators (members)
Matrix Matrix::operator+(Matrix M)
{
	return Plus(M);
}

Matrix Matrix::operator-(Matrix M)
{
	return Minus(M);
}

Matrix Matrix::operator*(const Matrix& M)
{
	return Product(M);
}

Vector Matrix::operator*(const Vector& v_rhs)
{
	return Product(v_rhs);
}

Matrix Matrix::operator*(double s)
{
	return Product(s);
}

Matrix Matrix::operator/(double s)
{
	return Division(s);
}

Matrix Matrix::operator=(Matrix M)
{
	components = M.components;
	rows	   = M.rows;
	columns	   = M.columns;
	return *this;
}

Matrix& Matrix::operator+=(const Matrix& M)
{
	if(rows != M.Columns() || columns != M.Rows())
	{
		std::cerr << "Error in Matrix::operator+=(): The dimensions of the two matrices to be added do not match." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		for(unsigned int i = 0; i < rows; i++)
			for(unsigned int j = 0; j < columns; j++)
				components[i][j] += M[i][j];
		return *this;
	}
}

Matrix& Matrix::operator-=(const Matrix& M)
{
	if(rows != M.Columns() || columns != M.Rows())
	{
		std::cerr << "Error in Matrix::operator-=(): The dimensions of the two matrices to be subtracted do not match." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		for(unsigned int i = 0; i < rows; i++)
			for(unsigned int j = 0; j < columns; j++)
				components[i][j] -= M[i][j];
		return *this;
	}
}

// Overloading operators (non-members)
Matrix operator*(double s, const Matrix& M)
{
	return M.Product(s);
}

Vector operator*(const Vector& v_left, const Matrix& M)
{
	if(v_left.Size() != M.Rows())
	{
		std::cerr << "Error in operator*(const Vector&,const Matrix&): The dimensions of the matrix and the vector do not match." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		std::vector<double> result_components(M.Columns(), 0.0);
		for(unsigned int i = 0; i < M.Columns(); i++)
			for(unsigned int j = 0; j < M.Rows(); j++)
				result_components[i] += v_left[j] * M[j][i];
		return Vector(result_components);
	}
}

std::ostream& operator<<(std::ostream& output, const Matrix& M)
{
	for(unsigned int i = 0; i < M.Rows(); i++)
	{
		if(i == 0)
			output << "\u2308";
		else if(i == M.Rows() - 1)
			output << "\u230A";
		else
			output << "|";
		for(unsigned int j = 0; j < M.Columns(); j++)
		{
			output << M[i][j];
			if(j < (M.Columns() - 1))
				output << "\t";
			else if(i == 0)
				output << "\u2309";
			else if(i == M.Rows() - 1)
				output << "\u230B";
			else
				output << "|";
		}
		if(i < M.Rows() - 1)
			output << std::endl;
	}
	return output;
}

bool operator==(const Matrix& M1, const Matrix& M2)
// Testing doubles for equality is still shitty of course.
{
	if(M1.Rows() != M2.Rows() || M1.Columns() != M2.Columns())
		return false;
	else
	{
		bool result = true;
		for(unsigned int i = 0; i < M1.Rows(); i++)
			for(unsigned int j = 0; j < M1.Columns(); j++)
				result *= (M1[i][j] == M2[i][j]);
		return result;
	}
}

Matrix Unit_Matrix(unsigned int dim)
{
	std::vector<double> ones(dim, 1.0);
	return Matrix(ones);
}

Matrix Rotation_Matrix(double alpha, int dim, Vector axis)
{
	double cosa = cos(alpha);
	double sina = sin(alpha);
	if(dim == 2)
	{
		std::vector<std::vector<double>> components =
			{
				{cosa, -sina},
				{sina, cosa}};
		return Matrix(components);
	}
	else if(dim == 3)
	{
		if(axis.Size() != 3)
		{
			std::cerr << "Error in Rotation_Matrix(): Rotation axis for 3D rotation is given via a " << dim << "D vector." << std::endl;
			std::exit(EXIT_FAILURE);
		}
		axis.Normalize();
		double n1 = axis[0];
		double n2 = axis[1];
		double n3 = axis[2];
		std::vector<std::vector<double>> components =
			{
				{cosa + n1 * n1 * (1.0 - cosa), n1 * n2 * (1.0 - cosa) - n3 * sina, n1 * n3 * (1.0 - cosa) + n2 * sina},
				{n1 * n2 * (1.0 - cosa) + n3 * sina, cosa + n2 * n2 * (1.0 - cosa), n2 * n3 * (1.0 - cosa) - n1 * sina},
				{n1 * n3 * (1.0 - cosa) - n2 * sina, n2 * n3 * (1.0 - cosa) + n1 * sina, cosa + n3 * n3 * (1.0 - cosa)}};
		return Matrix(components);
	}
	else
	{
		std::cerr << "Error in Rotation_Matrix(): Only 2D and 3D rotations are supported, not D = " << dim << "." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

}	// namespace libphysica
