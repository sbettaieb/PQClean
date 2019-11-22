#ifndef REPETITION_H
#define REPETITION_H

/**
 * @file repetition.h
 * @brief Header file for repetition.c
 */

#include <stdint.h>

void PQCLEAN_HQC1281CCA2_AVX2_LEAKTIME_repetition_code_encode(uint8_t *em, const uint8_t *m);
void PQCLEAN_HQC1281CCA2_AVX2_LEAKTIME_repetition_code_decode(uint8_t *m, const uint8_t *em);

#endif
