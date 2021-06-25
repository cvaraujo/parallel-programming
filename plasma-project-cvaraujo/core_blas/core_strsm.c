/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/core_blas/core_ztrsm.c, normal z -> s, Fri Sep 28 17:38:19 2018
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_trsm
 *
 *  Solves one of the matrix equations
 *
 *    \f[ op( A )\times X  = \alpha B, \f] or
 *    \f[ X \times op( A ) = \alpha B, \f]
 *
 *  where op( A ) is one of:
 *    \f[ op( A ) = A,   \f]
 *    \f[ op( A ) = A^T, \f]
 *    \f[ op( A ) = A^T, \f]
 *
 *  alpha is a scalar, X and B are m-by-n matrices, and
 *  A is a unit or non-unit, upper or lower triangular matrix.
 *  The matrix X overwrites B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          - PlasmaLeft:  op(A)*X = B,
 *          - PlasmaRight: X*op(A) = B.
 *
 * @param[in] uplo
 *          - PlasmaUpper: A is upper triangular,
 *          - PlasmaLower: A is lower triangular.
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          - PlasmaNonUnit: A has non-unit diagonal,
 *          - PlasmaUnit:    A has unit diagonal.
 *
 * @param[in] m
 *          The number of rows of the matrix B. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix B. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The lda-by-ka triangular matrix,
 *          where ka = m if side = PlasmaLeft,
 *            and ka = n if side = PlasmaRight.
 *          If uplo = PlasmaUpper, the leading k-by-k upper triangular part
 *          of the array A contains the upper triangular matrix, and the
 *          strictly lower triangular part of A is not referenced.
 *          If uplo = PlasmaLower, the leading k-by-k lower triangular part
 *          of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced.
 *          If diag = PlasmaUnit, the diagonal elements of A are also not
 *          referenced and are assumed to be 1.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[in,out] B
 *          On entry, the ldb-by-n right hand side matrix B.
 *          On exit, if return value = 0, the ldb-by-n solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void plasma_core_strsm(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                float alpha, const float *A, int lda,
                                                float *B, int ldb)
{
    cblas_strsm(CblasColMajor,
                (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)transa, (CBLAS_DIAG)diag,
                m, n,
                (alpha), A, lda,
                                    B, ldb);
}

/******************************************************************************/
void plasma_core_omp_strsm_default(
			   plasma_enum_t side, plasma_enum_t uplo,
			   plasma_enum_t transa, plasma_enum_t diag,
			   int m, int n,
			   float alpha, const float *A, int lda,
			   float *B, int ldb,
			   plasma_sequence_t *sequence, plasma_request_t *request)
{
  int ak;
  if (side == PlasmaLeft)
    ak = m;
  else
    ak = n;
  int side_ = side, uplo_ = uplo, transa_ = transa, diag_ = diag;
  int size_A = lda*ak, size_B = ldb*n;

#pragma omp task depend(in:A[0:lda*ak])		\
                     depend(inout:B[0:ldb*n])
  {
    plasma_core_strsm(side_, uplo_,
		      transa_, diag_,
		      m, n,
		      alpha, A, lda,
		      B, ldb);
  }
}

/******************************************************************************/
void plasma_core_omp_strsm(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n, float alpha,
    const float *A, int lda, float *B, int ldb,
    plasma_sequence_t *sequence, plasma_request_t *request)
{

  if (m <= 64 || n <= 64) {
    plasma_core_omp_strsm_default(side, uplo,
				  transa, diag,
				  m, n, alpha,
				  A, lda, B, ldb,
				  sequence, request);
  } else {
    int ak;
    if (side == PlasmaLeft) ak = m;
    else ak = n;

    int side_ = side, uplo_ = uplo, transa_ = transa, diag_ = diag;
    int size_A = lda * ak, size_B = ldb * n;
#pragma target nowait depend(in:A[0:lda*ak])				\
  depend(inout:B[0:ldb*n])						\
  map(tofrom: B[0:ldb*n])						\
  map(to: A[0:lda*ak])							\
  firstprivate(side_, uplo_, transa_, diag_, m, n, alpha, lda, ldb)

    if (sequence->status == PlasmaSuccess) {
      int transa_ = transa;
      int block_size = (lda < 128) ? lda : 128;
      int size_matrix = ldb;
      int m_new = m/block_size;
      int n_new = n/block_size;

#pragma omp parallel
#pragma omp single 
      {
	for(int m_ = 0; m_ < m_new; m_++) {
	  for(int n_ = 0; n_ < n_new; n_++) {
	    int k_ = m_;
	    int lda_ = m_ * block_size + k_ * block_size;
	    int ldb_ = k_ * block_size + n_ * block_size;
#pragma omp task
	    {
	      float A_new[block_size*block_size], B_new[block_size*block_size];
	      int i_a = m_, j_a = k_;
	      int i_b = k_, j_b = n_;
	      int insert = 0;

	      for(int l = 0; l < block_size; l++) {
		int bef_a = i_a * block_size + j_a * lda * block_size + l*size_matrix;
		int bef_b = i_b * block_size + j_b * lda * block_size + l*size_matrix;

		for (int k = 0; k < block_size; k++) {
		  A_new[insert] = A[bef_a++];
		  B_new[insert++] = B[bef_b++];
		}
	      }

	      plasma_core_strsm(side_, uplo_,
				transa_, diag_,
				block_size, block_size,
				alpha, (const float *) A_new, block_size,
				(float *) B_new, block_size);

	      insert = 0;
	      for(int l=0;l<block_size;l++) {
		int before_b = i_b*block_size+j_b*lda*block_size + l*size_matrix;
		for (int i = 0; i < block_size; i++) B[before_b++] = B_new[insert++];
	      }
	    }
	  }
	}
      }
    }
  }
}
