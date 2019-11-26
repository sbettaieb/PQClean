#ifndef VECTOR_H
#define VECTOR_H

/**
 * @file vector.h
 * @brief Header file for vector.c
 */

#include <stdint.h>
#include "nistseedexpander.h"

void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_set_random_fixed_weight_by_coordinates(AES_XOF_struct *ctx, uint32_t *v, const uint16_t weight);
void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_set_random_fixed_weight(AES_XOF_struct *ctx, uint8_t *v, const uint16_t weight);
void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_set_random(AES_XOF_struct *ctx, uint8_t *v);
void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_set_random_from_randombytes(uint8_t *v);

void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_add(uint8_t *o, const uint8_t *v1, const uint8_t *v2, const uint32_t size);
int PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_compare(const uint8_t *v1, const uint8_t *v2, const uint32_t size);

void PQCLEAN_HQC1921CCA2_REF_LEAKTIME_vect_resize(uint8_t *o, const uint32_t size_o, const uint8_t *v, const uint32_t size_v);

#endif
