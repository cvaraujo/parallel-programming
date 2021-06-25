/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/include/core_lapack_z.h, normal z -> c, Fri Sep 28 17:38:00 2018
 *
 **/

#ifndef PLASMA_CORE_LAPACK_C_H
#define PLASMA_CORE_LAPACK_C_H

#ifdef __cplusplus
extern "C" {
#endif

// LAPACK_GLOBAL is Fortran name mangling macro from LAPACKE

// LAPACKE_clantr broken (returns 0) in LAPACKE < 3.6.1
#ifndef LAPACK_clantr
#define LAPACK_clantr LAPACK_GLOBAL(clantr, CLANTR)
float LAPACK_clantr(const char *norm, const char *uplo, const char *diag,
                     const lapack_int *m, const lapack_int *n,
                     const plasma_complex32_t *A, const lapack_int *lda,
                     float *work);
#endif

// LAPACKE_clascl not available in LAPACKE < 3.6.0
#ifndef LAPACK_clascl
#define LAPACK_clascl LAPACK_GLOBAL(clascl, CLASCL)
void LAPACK_clascl(const char *type, const lapack_int *kl, const lapack_int *ku,
                   const float *cfrom, const float *cto,
                   const lapack_int *m, const lapack_int *n,
                   plasma_complex32_t *A, const lapack_int *lda,
                   lapack_int *info);
#endif

// LAPACKE_classq not available yet
#ifndef LAPACK_classq
#define LAPACK_classq LAPACK_GLOBAL(classq, CLASSQ)
void LAPACK_classq(const lapack_int *n, const plasma_complex32_t *x, const lapack_int *incx,
                   float *scale, float *sumsq);
#endif

// LAPACKE_clangb not available yet
#ifndef LAPACK_clangb
#define LAPACK_clangb LAPACK_GLOBAL(clangb, CLANGB)
float LAPACK_clangb(const char *norm,
                     const lapack_int *n, const lapack_int *kl, const lapack_int *ku,
                     const plasma_complex32_t *A, const lapack_int *lda,
                     float *work);

#endif

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // PLASMA_CORE_LAPACK_C_H
