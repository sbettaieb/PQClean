#ifndef BCH_H
#define BCH_H

/**
 * @file bch.h
 * Header file of bch.c
 */

#include <stddef.h>
#include <stdint.h>
#include "parameters.h"

void PQCLEAN_HQC2563CCA2_REF_LEAKTIME_bch_code_encode(uint8_t *codeword, const uint8_t *message);
void PQCLEAN_HQC2563CCA2_REF_LEAKTIME_bch_code_decode(uint8_t *message, uint8_t *vector);

#endif
