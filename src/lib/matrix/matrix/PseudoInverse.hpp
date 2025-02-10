/**
 * @file PseudoInverse.hpp
 *
 * Implementation of matrix pseudo inverse
 *
 * @author Julien Lecoeur <julien.lecoeur@gmail.com>
 * @author Julian Kent <julian@auterion.com>
 */

#pragma once

#include "math.hpp"

namespace matrix
{

/**
 * Geninv
 * Fast pseudoinverse based on full rank cholesky factorisation
 *
 * Courrieu, P. (2008). Fast Computation of Moore-Penrose Inverse Matrices, 8(2), 25–29. http://arxiv.org/abs/0804.4809
 */
template<typename Type, size_t M, size_t N>
bool geninv(const Matrix<Type, M, N> &G, Matrix<Type, N, M> &res)
{
	size_t rank;

	if (M <= N) {
		SquareMatrix<Type, M> A = G * G.transpose();
		SquareMatrix<Type, M> L = fullRankCholesky(A, rank);

		A = L.transpose() * L;
		SquareMatrix<Type, M> X;

		if (!inv(A, X, rank)) {
			res = Matrix<Type, N, M>();
			return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
		}

		// doing an intermediate assignment reduces stack usage
		A = X * X * L.transpose();

		auto temp_1 = L * A;
		auto temp_2 = G.transpose() * temp_1;
		res = temp_2;
		res = G.transpose() * (L * A);

	} else {
		SquareMatrix<Type, N> A = G.transpose() * G;
		SquareMatrix<Type, N> L = fullRankCholesky(A, rank);

		A = L.transpose() * L;
		SquareMatrix<Type, N> X;

		if (!inv(A, X, rank)) {
			res = Matrix<Type, N, M>();
			return false; // LCOV_EXCL_LINE -- this can only be hit from numerical issues
		}

		// doing an intermediate assignment reduces stack usage
		A = X * X * L.transpose();

		// Instead of eval()
		auto temp_1 = L * A;
		auto temp_2 = temp_1 * G.transpose();
		res = temp_2;
		// res = ((L * A) * G.transpose()).eval();
	}

	return true;
}


template<typename Type>
Type typeEpsilon();

template<> inline
float typeEpsilon<float>()
{
	return FLT_EPSILON;
}

/**
 * Full rank Cholesky factorization of A
 */
template<typename Type, size_t N>
SquareMatrix<Type, N> fullRankCholesky(const SquareMatrix<Type, N> &A,
				       size_t &rank)
{
	// Loses one ulp accuracy per row of diag, relative to largest magnitude
	const Type tol = N * typeEpsilon<Type>() * A.diag().max();

	// Initialize L to zero matrix
	Matrix<Type, N, N> L;
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			L(i, j) = 0;
		}
	}

	size_t r = 0;

	for (size_t k = 0; k < N; k++) {

		if (r == 0) {
			for (size_t i = k; i < N; i++) {
				L(i, r) = A(i, k);
			}

		} else {
			for (size_t i = k; i < N; i++) {
				// Compute LL = L[k:n, :r] * L[k, :r].T
				Type LL = Type();

				for (size_t j = 0; j < r; j++) {
					LL += L(i, j) * L(k, j);
				}

				L(i, r) = A(i, k) - LL;
			}
		}

		if (L(k, r) > tol) {
			L(k, r) = std::sqrt(L(k, r));

			if (k < N - 1) {
				for (size_t i = k + 1; i < N; i++) {
					L(i, r) = L(i, r) / L(k, r);
				}
			}

			r = r + 1;
		}
	}

	// Return rank
	rank = r;
	return L;
}

} // namespace matrix
