/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from /home/luszczek/workspace/plasma/bitbucket/plasma/core_blas/core_dcabs1.c, normal z -> c, Fri Sep 28 17:38:26 2018
 *
 **/

#include <plasma_core_blas.h>

#include <math.h>

/***************************************************************************//**
 *
 * @ingroup core_cabs1
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 *******************************************************************************
 *
 * @retval Complex 1-norm absolute value: abs(real(alpha)) + abs(imag(alpha)).
 *
 *******************************************************************************
 *
 * @sa plasma_core_scabs1
 *
 ******************************************************************************/
float plasma_core_scabs1(plasma_complex32_t alpha)
{
    return fabsf(creal(alpha)) + fabsf(cimag(alpha));
}
