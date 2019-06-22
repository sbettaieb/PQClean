#ifndef BCH_H
#define BCH_H

/**
 * \file bch.h
 * \brief Header file for bch.c
 */

#include <stdint.h>

#include "parameters.h"

void PQCLEAN_HQC1281CCA2_LEAKTIME_bch_code_encode(uint8_t *em, uint8_t *m);
void PQCLEAN_HQC1281CCA2_LEAKTIME_bch_code_decode(uint8_t *m, uint8_t *em);

#endif
