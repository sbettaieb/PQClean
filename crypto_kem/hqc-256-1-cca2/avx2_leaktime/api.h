#ifndef PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_API_H
#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_API_H

/**
 * \file api.h
 * \brief NIST KEM API used by the HQC_KEM IND-CCA2 scheme
 */

#include <stdint.h>

#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_CRYPTO_ALGNAME                      "HQC_256_1_CCA2"

#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_CRYPTO_SECRETKEYBYTES               8029
#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_CRYPTO_PUBLICKEYBYTES               7989
#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_CRYPTO_BYTES                        64
#define PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_CRYPTO_CIPHERTEXTBYTES              15961

// As a technicality, the public key is appended to the secret key in order to respect the NIST API.
// Without this constraint, CRYPTO_SECRETKEYBYTES would be defined as 32

int PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_crypto_kem_keypair(uint8_t *pk, uint8_t *sk);
int PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_crypto_kem_enc(uint8_t *ct, uint8_t *ss, const uint8_t *pk);
int PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_crypto_kem_dec(uint8_t *ss, const uint8_t *ct, const uint8_t *sk);

#endif
