/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/core_blas/core_zgemm.c, normal z -> s, Fri Sep 28 17:38:18 2018
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"

/***************************************************************************//**
 *
 * @ingroup core_gemm
 *
 *  Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^T, \f]
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m-by-k matrix, op( B ) a k-by-n matrix and C an m-by-n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transb
 *          - PlasmaNoTrans:   B is not transposed,
 *          - PlasmaTrans:     B is transposed,
 *          - PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transa = PlasmaNoTrans,
 *          and is m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          When transa = PlasmaNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transb = PlasmaNoTrans,
 *          and is k otherwise.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          When transb = PlasmaNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void plasma_core_sgemm(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
                float alpha, const float *A, int lda,
                                          const float *B, int ldb,
                float beta,        float *C, int ldc)
{
    cblas_sgemm(CblasColMajor,
                (CBLAS_TRANSPOSE)transa, (CBLAS_TRANSPOSE)transb,
                m, n, k,
                (alpha), A, lda,
                                    B, ldb,
                (beta),  C, ldc);
}

/******************************************************************************/
void plasma_core_omp_sgemm_default(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
    float alpha, const float *A, int lda,
                              const float *B, int ldb,
    float beta,        float *C, int ldc,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
  int ak;
  if (transa == PlasmaNoTrans)
    ak = k;
  else
    ak = m;

  int bk;
  if (transb == PlasmaNoTrans)
    bk = n;
  else
    bk = k;

#pragma omp task depend(in:A[0:lda*ak])		\
  depend(in:B[0:ldb*bk])			\
  depend(inout:C[0:ldc*n])
  {
    if (sequence->status == PlasmaSuccess)
      plasma_core_sgemm(transa, transb,
			m, n, k,
			alpha, A, lda,
			B, ldb,
			beta,  C, ldc);
  }
}

void plasma_core_omp_sgemm(
			   plasma_enum_t transa, plasma_enum_t transb,
			   int m, int n, int k,
			   float alpha, const float *A, int lda,
			   const float *B, int ldb,
			   float beta,        float *C, int ldc,
			   plasma_sequence_t *sequence, plasma_request_t *request) {

  
  if(m <= 64 || n <= 64) {
    plasma_core_omp_sgemm_default(transa, transb,
				  m, n, k,
				  alpha, A, lda,
				  B, ldb,
				  beta, C, ldc,
				  sequence, request);  
  } else {
    int ak;
    if (transa == PlasmaNoTrans)
      ak = k;
    else
      ak = m;

    int bk;
    if (transb == PlasmaNoTrans)
      bk = n;
    else
      bk = k;

    if (sequence->status == PlasmaSuccess) {
      int transa_ = transa, transb_ = transb;
      int size_A = lda*ak, size_B = ldb*bk,size_C = ldc*n;

#pragma omp target nowait			\
  depend(in:A[0:lda*ak])			\
  depend(in:B[0:ldb*bk])			\
  depend(inout:C[0:ldc*n])			\
  firstprivate(m,n,k,alpha,beta,lda,ak)		\
  firstprivate(ldb,bk,ldc,transa_,transb_)	\
  map(to:A[0:size_A],B[0:size_B])		\
  map(tofrom:C[0:size_C])
      {

	int block_size = (lda < 128) ? lda : 128;
	int size_matrix = lda;
	int m_new = m/block_size;
	int n_new = n/block_size;
	int k_new = k/block_size;
	int zbeta;

#pragma omp parallel
#pragma omp single
	{
	  for (int m_ = 0; m_ < m_new; m_++){
	    for (int n_ = 0; n_ < n_new; n_++) {
	      for (int k_ = 0; k_ < k_new; k_++) {
		int lda_ = m_ * block_size + k_ * block_size;
		int ldb_ = k_ * block_size + n_ * block_size;
		int ldc_ = m_ * block_size + n_ * block_size;
#pragma omp task							\
  depend(in:A[0:lda_])							\
  depend(in:B[0:ldb_])							\
  depend(inout:C[m_*block_size+n_*lda*block_size:m_*block_size+n_*lda*block_size + block_size * block_size])
		{

		  float A_new[block_size*block_size], B_new[block_size*block_size], C_new[block_size*block_size];
		  int i_a = m_, j_a=k_;
		  int i_b = k_, j_b=n_;
		  int i_c = m_, j_c=n_;
		  int insert = 0;

		  for(int l = 0; l < block_size; l++) {
		    int before_a = i_a*block_size+j_a*lda*block_size + l*size_matrix;
		    int before_b = i_b*block_size+j_b*lda*block_size + l*size_matrix;
		    int before_c = i_c*block_size+j_c*lda*block_size + l*size_matrix;

		    for (int i = 0; i < block_size; i++) {
		      A_new[insert] = A[before_a++];
		      B_new[insert] = B[before_b++];
		      C_new[insert++] = C[before_c++];
		    }
		  }

		  float zbeta = k_== 0 ? beta : 1.0;

		  plasma_core_sgemm(transa_, transb_,
				    block_size, block_size, block_size,
				    alpha, (const float *) A_new, block_size,
				    (const float *) B_new, block_size,
				    zbeta, (float *) C_new, block_size);

		  insert = 0;
		  for(int l = 0; l < block_size; l++) {
		    int before_c = i_c*block_size+j_c*lda*block_size + l*size_matrix;

		    for (int i = 0; i < block_size; i++) C[before_c++]=C_new[insert++];

		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
