#include <iostream>
#include <cmath>
using namespace std;

//10.15.1 C++ code: LU decomposition
int LU_decomposition(const int n,
					std::vector<std::vector<double>> & a,
					std::vector<int> & swap_indices,
					int & swap_count)
{
	const double tol = 1.0e-14; // tolerance for pivot

	swap_indices.clear();
	swap_count = 0;

	if (n < 1) return 1; // fail
	if (a.size() < n) return 1; // fail

	int i=0;
	int j=0;

	swap_indices.reserve(n);
	for (i = 0; i < n; ++i) {
		swap_indices.push_back(i); // initialize to (0,...n,-1)
	}

	for (i = 0; i < n; ++i) {
	// step 1: calculate scaled coeffs, compare for max pivot
		int pivot_row = i;
		double max_pivot = 0;
		for (j=i; j < n; ++j) {
			int k = i;
			double a_hat = std::abs(a[j][k]);
			if (a_hat > 0.0) {
				for (k=i+1; k < n; ++k) {
					double tmp = std::abs(a[j][k]);
					if (a_hat < tmp) a_hat = tmp;
				}
				double a_scaled = std::abs(a[j][i])/a_hat; // note: a_ji not a_ij
				if (max_pivot < a_scaled) {
					max_pivot = a_scaled;
					pivot_row = j;
				}
			}
		}

		// step 2: if max pivot <= tol return fail (inconsistent or not linearly independent)
		if (std::abs(a[pivot_row][i]) <= tol) {
			swap_indices.clear();
			swap_count = 0;
			return 1; // fail
		}
		// step 3: found the max pivot, swap rows if row != i
		if (pivot_row != i) {
			// swap the whole row, including lower triangular entries
			// keep track of number of swaps
			// update the swap index
			++swap_count;
			int si = swap_indices[i];
			swap_indices[i] = swap_indices[pivot_row];
			swap_indices[pivot_row] = si;

			for (j=0; j < n; ++j) {
				double tmp = a[i][j];
				a[i][j] = a[pivot_row][j];
				a[pivot_row][j] = tmp;
			}
		}
		double inv_pivot = 1.0/a[i][i]; // this is nonzero

		// step 4: eliminate column i from rows i+1 <= j < n, overwrite make LU matrix
		for (j=i+1; j < n; ++j) {
			double multiplier = a[j][i]*inv_pivot;
			a[j][i] = multiplier; // store multiplier in lower triangular matrix
	
			for (int k=i+1; k < n; ++k) {
				a[j][k] -= multiplier*a[i][k]; // calculate upper triangular matrix
			}
		}
	}
	return 0;
}

//10.15.2 C++ code: LU solver
int LU_solve(const int n,
			const std::vector<std::vector<double>> & LU,
			const std::vector<int> & swap_indices,
			const int swap_count,
			const std::vector<double> & rhs,
			std::vector<double> & x)
	{
		x.clear();
		if (n < 1) return 1; // fail

		std::vector<double> y(n, 0.0); // temporary storage
		// forward substitution Ly = rhs
		int i=0;
		for (i = 0; i < n; ++i) {
			int si = swap_indices[i];
			double sum = rhs[si];
			for (int j = 0; j < i; ++j) {
				sum -= LU[i][j] * y[j];
			}
		y[i] = sum;
	}
	// backward substitution Ux = y
	for (i = n-1; i >= 0; --i) {
		double sum = y[i];
		for (int j = i+1; j < n; ++j) {
			sum -= LU[i][j] * x[j];
		}
		x[i] = sum / LU[i][i];
	}
	return 0;
}
	


//10.15.3 C++ code: LU inverse
int LU_inverse(const int n,
				const std::vector<std::vector<double>> & LU,
				const std::vector<int> & swap_indices,
				const int swap_count,
				std::vector<std::vector<double>> & inv_matrix)
{
	if (n < 1) return 1; // fail
	for (int j = 0; j < n; ++j) {
		std::vector<double> rhs(n, 0.0);
		std::vector<double> x(n, 0.0);
		rhs[j] = 1.0;

		inv_matrix[j].clear();

		int rc = LU_solve(n, LU, swap_indices, swap_count, rhs, x);
		if (rc) {
			for (j = 0; j < n; ++j) {
				inv_matrix[j].clear(); // fail, clear everything
			}
			return rc; // fail
		}
		for (int i = 0; i < n; ++i) {
			inv_matrix[i][j] = x[i];
		}
	}
	return 0;
}

//10.15.4 C++ code: LU determinant
int LU_determinant(const int n,
					const std::vector<std::vector<double>> & LU,
					const std::vector<int> & swap_indices,
					const int swap_count,
					double & det)
{
	det = 0;
	if (n < 1) return 1; // fail
	det = ((swap_count % 2) == 0) ? 1 : -1;
	for (int i = 0; i < n; ++i) {
		det *= LU[i][i];
	}
	return 0;
}