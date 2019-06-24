#ifndef BCH_H
#define BCH_H

/**
 * \file bch.h
 * \brief Header file for bch.c
 */

#include "parameters.h"

#include <stdint.h>

void PQCLEAN_HQC1281CCA2_LEAKTIME_bch_code_encode(uint8_t *em, const uint8_t *m);
void PQCLEAN_HQC1281CCA2_LEAKTIME_bch_code_decode(uint8_t *m, const uint8_t *em);

#endif
