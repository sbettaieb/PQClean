/**
 * \file kem.c
 * \brief Implementation of api.h
 */

#include <stdint.h>
#include <string.h>

#include "api.h"

#include "hqc.h"
#include "parameters.h"
#include "parsing.h"
#include "vector.h"

#include "nistseedexpander.h"
#include "sha2.h"

/**
 *\fn int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_keypair(unsigned char* pk, unsigned char* sk)
 *\brief Keygen of the HQC_KEM IND_CAA2 scheme
 *
 * The public key is composed of the syndrome <b>s</b> as well as the seed used to generate the vector <b>h</b>.
 *
 * The secret key is composed of the seed used to generate vectors <b>x</b> and <b>y</b>.
 * As a technicality, the public key is appended to the secret key in order to respect NIST API.
 *
 * \param[out] pk String containing the public key
 * \param[out] sk String containing the secret key
 * return 0 if keygen is successful
 */
int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_keypair(unsigned char *pk, unsigned char *sk) {
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_pke_keygen(pk, sk);
    return 0;
}

/**
 *\fn int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk)
 *\brief Encapsulation of the HQC_KEM IND_CAA2 scheme
 *
 * \param[out] ct String containing the ciphertext
 * \param[out] ss String containing the shared secret
 * \param[in] pk String containing the public key
 * return 0 if encapsulation is successful
 */
int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk) {
    // Computing m
    uint8_t m[VEC_K_SIZE_BYTES] = {0};
    PQCLEAN_HQC1281CCA2_LEAKTIME_vect_set_random_from_randombytes(m);

    // Generating G function
    unsigned char diversifier_bytes[8];
    for (int i = 0; i < 8; ++i) {
        diversifier_bytes[i] = 0;
    }
    unsigned char seed_G[VEC_K_SIZE_BYTES];
    memcpy(seed_G, m, VEC_K_SIZE_BYTES);
    AES_XOF_struct G_seedexpander;
    seedexpander_init(&G_seedexpander, seed_G, diversifier_bytes, SEEDEXPANDER_MAX_LENGTH);

    // Computing theta
    unsigned char theta[SEED_BYTES];
    seedexpander(&G_seedexpander, theta, SEED_BYTES);

    // Encrypting m
    uint8_t u [VEC_N_SIZE_BYTES] = {0};
    uint8_t v [VEC_N1N2_SIZE_BYTES] = {0};
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_pke_encrypt(u, v, m, theta, pk);

    // Computing d
    unsigned char d[SHA512_BYTES];
    sha512(d, m, VEC_K_SIZE_BYTES);

    // Computing shared secret
    unsigned char mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES];
    memcpy(mc, m, VEC_K_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES, u, VEC_N_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES, v, VEC_N1N2_SIZE_BYTES);
    sha512(ss, mc, VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES);

    // Computing ciphertext
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_ciphertext_to_string(ct, u, v, d);

    return 0;
}

/**
 *\fn int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk)
 *\brief Decapsulation of the HQC_KEM IND_CAA2 scheme
 *
 * \param[out] ss String containing the shared secret
 * \param[in] ct String containing the cipĥertext
 * \param[in] sk String containing the secret key
 * return 0 if decapsulation is successful, -1 otherwise
 */
int PQCLEAN_HQC1281CCA2_LEAKTIME_crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk) {
    // Retrieving u, v and d from ciphertext
    uint8_t u [VEC_N_SIZE_BYTES] = {0};
    uint8_t v [VEC_N1N2_SIZE_BYTES] = {0};
    unsigned char d[SHA512_BYTES];
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_ciphertext_from_string(u, v, d, ct);

    // Retrieving pk from sk
    unsigned char pk[PUBLIC_KEY_BYTES];
    memcpy(pk, sk + SEED_BYTES, PUBLIC_KEY_BYTES);

    //Decryting
    uint8_t m[VEC_K_SIZE_BYTES] = {0};
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_pke_decrypt(m, u, v, sk);

    // Generating G function
    unsigned char diversifier_bytes[8];
    for (int i = 0; i < 8; ++i) {
        diversifier_bytes[i] = 0;
    }
    unsigned char seed_G[VEC_K_SIZE_BYTES];
    memcpy(seed_G, m, VEC_K_SIZE_BYTES);
    AES_XOF_struct G_seedexpander;
    seedexpander_init(&G_seedexpander, seed_G, diversifier_bytes, SEEDEXPANDER_MAX_LENGTH);

    // Computing theta
    unsigned char theta[SEED_BYTES];
    seedexpander(&G_seedexpander, theta, SEED_BYTES);

    // Encrypting m'
    uint8_t u2 [VEC_N_SIZE_BYTES] = {0};
    uint8_t v2 [VEC_N1N2_SIZE_BYTES] = {0};
    PQCLEAN_HQC1281CCA2_LEAKTIME_hqc_pke_encrypt(u2, v2, m, theta, pk);

    // Checking that c = c' and abort otherwise
    int abort = 0;
    if (PQCLEAN_HQC1281CCA2_LEAKTIME_vect_compare(u, u2, VEC_N_SIZE_BYTES) != 0) {
        abort = 1;
    }
    if (PQCLEAN_HQC1281CCA2_LEAKTIME_vect_compare(v, v2, VEC_N1N2_SIZE_BYTES) != 0) {
        abort = 1;
    }

    // Computing d'
    unsigned char d2[SHA512_BYTES];
    sha512(d2, m, VEC_K_SIZE_BYTES);

    // Checking that d = d' and abort otherwise
    if (memcmp(d, d2, SHA512_BYTES) != 0) {
        abort = 1;
    }

    if (abort == 1) {
        memset(ss, 0, SHARED_SECRET_BYTES);
        return -1;
    }

    // Computing shared secret
    unsigned char mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES];
    memcpy(mc, m, VEC_K_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES, u, VEC_N_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES, v, VEC_N1N2_SIZE_BYTES);
    sha512(ss, mc, VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES);

    return 0;
}
