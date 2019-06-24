#ifndef PARSING_H
#define PARSING_H

/**
 * \file parsing.h
 * \brief Header file for parsing.c
 */

#include <stdint.h>

void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_secret_key_to_string(unsigned char *sk, const unsigned char *sk_seed, const unsigned char *pk);
void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_secret_key_from_string(uint8_t *x, uint8_t *y, unsigned char *pk, const unsigned char *sk);

void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_public_key_to_string(unsigned char *pk, const unsigned char *pk_seed, const uint8_t *s);
void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_public_key_from_string(uint8_t *h, uint8_t *s, const unsigned char *pk);

void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_ciphertext_to_string(unsigned char *ct, const uint8_t *u, const uint8_t *v, const unsigned char *d);
void PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_ciphertext_from_string(uint8_t *u, uint8_t *v, unsigned char *d, const unsigned char *ct);

#endif
