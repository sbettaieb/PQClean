#ifndef HQC_PARAMETERS_H
#define HQC_PARAMETERS_H

/**
 * \file parameters.h
 * \brief Parameters of the HQC_KEM IND-CCA2 scheme
 */

#include "api.h"

#define CEIL_DIVIDE(a, b)  (((a)/(b)) + ((a) % (b) == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/

/*
  #define PARAM_N                             Define the parameter n of the scheme
  #define PARAM_N1                            Define the parameter n1 of the scheme (length of BCH code)
  #define PARAM_N2                            Define the parameter n2 of the scheme (length of the repetition code)
  #define PARAM_N1N2                          Define the parameter n1 * n2 of the scheme (length of the tensor code)
  #define PARAM_T                             Define a threshold for decoding repetition code word (PARAM_T = (PARAM_N2 - 1) / 2)
  #define PARAM_OMEGA                         Define the parameter omega of the scheme
  #define PARAM_OMEGA_E                       Define the parameter omega_e of the scheme
  #define PARAM_OMEGA_R                       Define the parameter omega_r of the scheme
  #define PARAM_DELTA                         Define the parameter delta of the scheme (correcting capacity of the BCH code)
  #define PARAM_SECURITY                      Define the security level corresponding to the chosen parameters
  #define PARAM_DFR_EXP                       Define the decryption failure rate corresponding to the chosen parameters

  #define SECRET_KEY_BYTES                    Define the size of the secret key in bytes
  #define PUBLIC_KEY_BYTES                    Define the size of the public key in bytes
  #define SHARED_SECRET_BYTES                 Define the size of the shared secret in bytes
  #define CIPHERTEXT_BYTES                    Define the size of the ciphertext in bytes

  #define PARAM_POLY                          Define p(x) = 1 + x^3 + x^10 the primitive polynomial of degree 10 in hexadecimal representation
                                              We define alpha as the root of p(x), i.e p(alpha) = 0. Then, alpha is a primitive element,
                                              and the powers of alpha generate all the nonzero elements of GF(2^10)
                                              We use this polynomial to build the Galois Field GF(2^10) needed for the BCH code
  #define PARAM_M                             Define a positive integer
  #define PARAM_GF_MUL_ORDER                  Define the size of the multiplicative group of GF(2^m), i.e 2^m -1
  #define PARAM_G                             Define the size of the generator polynomial of BCH code
  #define PARAM_K                             Define the size of the information bits of the BCH code

  #define UTILS_REJECTION_THRESHOLD           Define the rejection threshold used to generate given weight vectors (see vector_u32_fixed_weight function documentation)

  #define VEC_N_ARRAY_SIZE_BYTES              Define the size of the array used to store a PARAM_N sized vector in bytes
  #define VEC_N1_ARRAY_SIZE_BYTES             Define the size of the array used to store a PARAM_N1 sized vector in bytes
  #define VEC_N1N2_ARRAY_SIZE_BYTES           Define the size of the array used to store a PARAM_N1N2 sized vector in bytes
  #define VEC_K_ARRAY_SIZE_BYTES              Define the size of the array used to store a PARAM_K sized vector in bytes

  #define SHA512_BYTES                        Define the size of SHA512 output in bytes
  #define SEED_BYTES                          Define the size of the seed in bytes
  #define SEEDEXPANDER_MAX_LENGTH             Define the seed expander max length
*/

#define PARAM_N                             43669
#define PARAM_N1                            766
#define PARAM_N2                            57
#define PARAM_N1N2                          43662
#define PARAM_T                             28
#define PARAM_OMEGA                         101
#define PARAM_OMEGA_E                       117
#define PARAM_OMEGA_R                       117
#define PARAM_DELTA                         57
#define PARAM_SECURITY                      192
#define PARAM_DFR_EXP                       128

#define SECRET_KEY_BYTES                    PQCLEAN_HQC1921CCA2_LEAKTIME_CRYPTO_SECRETKEYBYTES
#define PUBLIC_KEY_BYTES                    PQCLEAN_HQC1921CCA2_LEAKTIME_CRYPTO_PUBLICKEYBYTES
#define SHARED_SECRET_BYTES                 PQCLEAN_HQC1921CCA2_LEAKTIME_CRYPTO_BYTES
#define CIPHERTEXT_BYTES                    PQCLEAN_HQC1921CCA2_LEAKTIME_CRYPTO_CIPHERTEXTBYTES

#define PARAM_POLY                          0x409
#define PARAM_M                             10
#define PARAM_GF_MUL_ORDER                  1023
#define PARAM_G                             511
#define PARAM_K                             256

#define UTILS_REJECTION_THRESHOLD           16768896

#define VEC_N_SIZE_BYTES                    CEIL_DIVIDE(PARAM_N, 8)
#define VEC_N1_SIZE_BYTES                   CEIL_DIVIDE(PARAM_N1, 8)
#define VEC_N1N2_SIZE_BYTES                 CEIL_DIVIDE(PARAM_N1N2, 8)
#define VEC_K_SIZE_BYTES                    CEIL_DIVIDE(PARAM_K, 8)

#define SHA512_BYTES                        64
#define SEED_BYTES                          40
#define SEEDEXPANDER_MAX_LENGTH             4294967295

#define BCH_GEN_POLY { \
        0xC2, 0x6D, 0x64, 0xFC, 0xDE, 0xF4, 0x4E, 0xD7, 0x52, 0x1E, 0xF8, 0x49, \
        0xC6, 0x51, 0x08, 0xB0, 0xCD, 0x57, 0xA6, 0xB3, 0x7D, 0x74, 0x6F, 0x78, \
        0x40, 0x83, 0xCD, 0x09, 0x29, 0xAF, 0xB5, 0xF8, 0xBF, 0x56, 0x33, 0x7F, \
        0x1E, 0x89, 0xB2, 0xD5, 0x83, 0xFF, 0x6B, 0xFB, 0x62, 0x7D, 0x42, 0xF4, \
        0x12, 0xF8, 0x1F, 0x49, 0xA0, 0x18, 0xE6, 0x34, 0x8A, 0x8F, 0x89, 0xB2, \
        0xDC, 0x37, 0x42, 0x26 \
    };
#define BCH_MSG_MASK_M1                     0xC0
#define BCH_MSG_MASK_M2                     0x3F

#endif
