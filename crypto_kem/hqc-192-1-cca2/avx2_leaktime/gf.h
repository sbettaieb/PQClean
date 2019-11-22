#ifndef GF_H
#define GF_H

/**
 * @file gf.h
 * Header file of gf.c
 */

#include <stddef.h>
#include <stdint.h>

uint16_t PQCLEAN_HQC1921CCA2_AVX2_LEAKTIME_gf_log(const uint16_t elt);
uint16_t PQCLEAN_HQC1921CCA2_AVX2_LEAKTIME_gf_mul(const uint16_t a, const uint16_t b);
uint16_t PQCLEAN_HQC1921CCA2_AVX2_LEAKTIME_gf_square(const uint16_t a);
uint16_t PQCLEAN_HQC1921CCA2_AVX2_LEAKTIME_gf_inverse(const uint16_t a);
uint16_t PQCLEAN_HQC1921CCA2_AVX2_LEAKTIME_gf_mod(const uint16_t i);

#endif
