// matrix_test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

class SplineInterpolator
{
private:
	gsl_vector* input_data, * rhs, * D;
	gsl_vector* a, * b, * c, * d;	// Spline parameters
	gsl_matrix* lhs;

	void store_input_data(const double* data, int data_size);
	void initialize_lhs();
	void initialize_rhs();
	void solve_for_D();
	void calculate_spline_parameters();

public:
	SplineInterpolator(const double* data, int data_size);
	void write_lhs_data(std::string file_name);
	void write_rhs_data(std::string file_name);
	void write_D_data(std::string file_name);

	enum BoundaryCondition
	{
		Clamped,
		Natural,
		Periodic
	};
};

int main()
{
    std::cout << "Hello World!\n";

	const std::vector<double> data{30.0, 45.0, 15.0, 22.0, 30.0, 68.0, 133.0, 90.0, 59.0, 38.0, 19.0, 10.0, 30.0};

	SplineInterpolator spline(data.data(), data.size());
	spline.write_lhs_data("lhs_matrix.txt");
	spline.write_rhs_data("rhs_vector.txt");
	spline.write_D_data("D_vector.txt");
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

SplineInterpolator::SplineInterpolator(const double* data, int data_size)
{
	// Load input data
	store_input_data(data, data_size);

	// Initialize the left-hand-side matrix
	initialize_lhs();

	// Calculate the right-hand-side vector
	initialize_rhs();

	// Solve for D values
	solve_for_D();
}

void SplineInterpolator::store_input_data(const double* data, int data_size)
{
	// Initialize new vector to store data
	input_data = gsl_vector_calloc(data_size);

	// Load input data
	for (int i = 0; i < data_size; i++)
	{
		gsl_vector_set(input_data, i, data[i]);
	}
}

void SplineInterpolator::initialize_lhs()
{
	// Get data size
	int size = input_data->size;

	// Initialize matrix with 0 values
	lhs = gsl_matrix_calloc(size, size);

	// Set values for the diagonal
	gsl_matrix_add_diagonal(lhs, 4.0);

	// Set values for the (0,0) and (size-1,size-1) positions
	gsl_matrix_set(lhs, 0, 0, 2.0);
	gsl_matrix_set(lhs, size - 1, size - 1, 2.0);

	// Set values for the super-diagonal
	gsl_vector_view m1_off_diag = gsl_matrix_superdiagonal(lhs, 1);
	gsl_vector_set_all(&m1_off_diag.vector, 1.0);

	// Set values for the sub-diagonal
	m1_off_diag = gsl_matrix_subdiagonal(lhs, 1);
	gsl_vector_set_all(&m1_off_diag.vector, 1.0);
}

void SplineInterpolator::initialize_rhs()
{
	// Get data size
	int size = input_data->size;

	// Initialize right-hand-side vector with zero values
	rhs = gsl_vector_calloc(size);

	// Calculate values for right-hand-side vector
	double val1, val2, current_value;
	for (int i = 0; i < size; i++)
	{
		// Calculate intermediate values
		val1 = gsl_vector_get(input_data, std::min(i + 1, size - 1));
		val2 = gsl_vector_get(input_data, std::max(i - 1, 0));
		current_value = 3.0 * (val1 - val2);

		// Store final value
		gsl_vector_set(rhs, i, current_value);
	}
}

void SplineInterpolator::solve_for_D()
{
	// Get data size
	int size = input_data->size;

	// Create LU and initialize to values in left-hand-side matrix
	gsl_matrix* LU = gsl_matrix_alloc(size, size);
	gsl_matrix_memcpy(LU, lhs);

	// Create permutation vector to store permutation
	gsl_permutation* p = gsl_permutation_calloc(size);

	// Create pointer to int for third argument to LU decomposition function
	int* signum = new int(0);

	// Get the LU decomposition
	gsl_linalg_LU_decomp(LU, p, signum);

	// Initialize D (to store results of LU solve)
	D = gsl_vector_calloc(size);

	// Solve for D
	const gsl_vector* rhs_const = rhs;
	gsl_linalg_LU_solve(LU, p, rhs, D);
}

void SplineInterpolator::calculate_spline_parameters()
{
	// Store spline parameters 'a'
	gsl_vector_memcpy(a, input_data);

	// Store spline parameters 'b'
	gsl_vector_memcpy(b, D);


}

void SplineInterpolator::write_lhs_data(std::string file_name)
{
	// Write contents of matrix to a file
	FILE* file;
	const gsl_matrix* lhs_const = lhs;
	errno_t err = fopen_s(&file, file_name.data(), "w");
	gsl_matrix_fprintf(file, lhs_const, "%f");
	if (file) fclose(file);
}

void SplineInterpolator::write_rhs_data(std::string file_name)
{
	// Write contents of matrix to a file
	FILE* file;
	const gsl_vector* rhs_const = rhs;
	errno_t err = fopen_s(&file, file_name.data(), "w");
	gsl_vector_fprintf(file, rhs_const, "%f");
	if (file) fclose(file);
}

void SplineInterpolator::write_D_data(std::string file_name)
{
	// Write contents of matrix to a file
	FILE* file;
	const gsl_vector* D_const = D;
	errno_t err = fopen_s(&file, file_name.data(), "w");
	gsl_vector_fprintf(file, D_const, "%f");
	if (file) fclose(file);
}

