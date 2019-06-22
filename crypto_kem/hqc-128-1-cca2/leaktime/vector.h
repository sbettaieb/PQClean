#ifndef VECTOR_H
#define VECTOR_H

/**
 * \file vector.h
 * \brief Header file for vector.c
 */

#include <stdint.h>

#include "nistseedexpander.h"

void PQCLEAN_HQC1281CCA2_LEAKTIME_vect_fixed_weight(uint8_t *v, const uint16_t weight, AES_XOF_struct *ctx);
void PQCLEAN_HQC1281CCA2_LEAKTIME_vect_set_random(uint8_t *v, AES_XOF_struct *ctx);
void PQCLEAN_HQC1281CCA2_LEAKTIME_vect_set_random_from_randombytes(uint8_t *v);

void PQCLEAN_HQC1281CCA2_LEAKTIME_vect_add(uint8_t *o, uint8_t *v1, uint8_t *v2, uint32_t size);
int PQCLEAN_HQC1281CCA2_LEAKTIME_vect_compare(uint8_t *v1, uint8_t *v2, uint32_t size);

void PQCLEAN_HQC1281CCA2_LEAKTIME_vect_resize(uint8_t *o, uint32_t size_o, uint8_t *v, uint32_t size_v);

#endif
