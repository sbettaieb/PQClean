#ifndef GF2X_H
#define GF2X_H

/**
 * @file gf2x.h
 * @brief Header file for gf2x.c
 */

#include <stdint.h>

void PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_vect_mul(uint8_t *o, const uint32_t *v1, const uint8_t *v2, const uint32_t weight);

#endif
