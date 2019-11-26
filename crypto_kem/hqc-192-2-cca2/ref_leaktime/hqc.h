#ifndef HQC_H
#define HQC_H

/**
 * @file hqc.h
 * @brief Functions of the HQC_PKE IND_CPA scheme
 */

#include <stdint.h>

void PQCLEAN_HQC1922CCA2_REF_LEAKTIME_hqc_pke_keygen(uint8_t *pk, uint8_t *sk);
void PQCLEAN_HQC1922CCA2_REF_LEAKTIME_hqc_pke_encrypt(uint8_t *u, uint8_t *v, const uint8_t *m, const uint8_t *theta, const uint8_t *pk);
void PQCLEAN_HQC1922CCA2_REF_LEAKTIME_hqc_pke_decrypt(uint8_t *m, const uint8_t *u, const uint8_t *v, const uint8_t *sk);

#endif
