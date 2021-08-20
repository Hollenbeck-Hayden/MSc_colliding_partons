#ifndef MATRIX_H
#define MATRIX_H

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

using matrix_type = std::vector<std::vector<double>>;

void print_matrix(const matrix_type& m)
{
	for (const auto& row : m)
	{
		for (const double x : row)
			std::cout << std::setw(8) << std::setprecision(2) << (fabs(x) < 1e-24 ? 0 : x) << "  ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

matrix_type mat_mult(const matrix_type& a, const matrix_type& b)
{
	if (a[0].size() != b.size()) {
		std::cout << "ERROR: WRONG MATRIX MULT SIZE - " << a.size() << "x" << a[0].size() << " times " << b.size() << "x" << b[0].size() << std::endl;
	}
	matrix_type c(a.size(), std::vector<double>(b[0].size(), 0));
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < b[0].size(); j++)
			for (int k = 0; k < a[0].size(); k++)
				c[i][j] += a[i][k] * b[k][j];
	return c;
}

matrix_type transpose(const matrix_type& a)
{
	matrix_type t(a.size(), std::vector<double>(a[0].size(), 0));
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[0].size(); j++)
			t[j][i] = a[i][j];
	return t;
}

/*
 * Cholesky LDLT inverse of a real, symmetric, positive-definite
 * matrix (equations from wikipedia).
 */
std::vector<std::vector<double>> cholesky_inverse(const std::vector<std::vector<double>>& A)
{
	const auto dim = A.size();
	std::vector<std::vector<double>> L(dim, std::vector<double>(dim, 0));
	std::vector<std::vector<double>> D(dim, std::vector<double>(dim, 0));

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double total = 0;
			for (int k = 0; k < j; k++)
				total += L[i][k] * L[j][k] * D[k][k];

			if (i == j)
			{
				D[i][j] = A[i][j] - total;
				L[i][j] = 1;
			}
			else
			{
				L[i][j] = (A[i][j] - total) / D[j][j];
			}
		}
	}

	matrix_type L_inv(dim, std::vector<double>(dim, 0));
	matrix_type D_inv(dim, std::vector<double>(dim, 0));
	matrix_type N(dim, std::vector<double>(dim, 0));

	for (int i = 0; i < dim; i++)
		for (int j = 0; j < i; j++)
			N[i][j] = -L[i][j];

	const matrix_type N_original = N;

	for (int i = 0; i < dim; i++)
	{
		L_inv[i][i] = 1;
		D_inv[i][i] = 1./D[i][i];
	}
	
	for (int k = 0; k < dim; k++)
	{
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < i; j++)
					L_inv[i][j] += N[i][j];
		N = mat_mult(N, N_original);
	}

	return mat_mult(transpose(L_inv), mat_mult(D_inv, L_inv));
}

#endif
