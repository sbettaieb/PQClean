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

#define PARAM_N                             70853
#define PARAM_N1                            796
#define PARAM_N2                            89
#define PARAM_N1N2                          70844
#define PARAM_T                             44
#define PARAM_OMEGA                         133
#define PARAM_OMEGA_E                       153
#define PARAM_OMEGA_R                       153
#define PARAM_DELTA                         60
#define PARAM_SECURITY                      256
#define PARAM_DFR_EXP                       256

#define SECRET_KEY_BYTES                    PQCLEAN_HQC2563CCA2_LEAKTIME_CRYPTO_SECRETKEYBYTES
#define PUBLIC_KEY_BYTES                    PQCLEAN_HQC2563CCA2_LEAKTIME_CRYPTO_PUBLICKEYBYTES
#define SHARED_SECRET_BYTES                 PQCLEAN_HQC2563CCA2_LEAKTIME_CRYPTO_BYTES
#define CIPHERTEXT_BYTES                    PQCLEAN_HQC2563CCA2_LEAKTIME_CRYPTO_CIPHERTEXTBYTES

#define PARAM_POLY                          0x409
#define PARAM_M                             10
#define PARAM_GF_MUL_ORDER                  1023
#define PARAM_G                             541
#define PARAM_K                             256


#define UTILS_REJECTION_THRESHOLD           16721308

#define VEC_N_SIZE_BYTES                    CEIL_DIVIDE(PARAM_N, 8)
#define VEC_N1_SIZE_BYTES                   CEIL_DIVIDE(PARAM_N1, 8)
#define VEC_N1N2_SIZE_BYTES                 CEIL_DIVIDE(PARAM_N1N2, 8)
#define VEC_K_SIZE_BYTES                    CEIL_DIVIDE(PARAM_K, 8)

#define SHA512_BYTES                        64
#define SEED_BYTES                          40
#define SEEDEXPANDER_MAX_LENGTH             4294967295

#define BCH_GEN_POLY { \
        0xDA, 0x15, 0xFB, 0x7C, 0x96, 0x9C, 0xE2, 0xAC, 0xA4, 0x7F, 0x3D, 0x5A, \
        0x5C, 0x78, 0xBD, 0x3D, 0x13, 0xE8, 0xB1, 0x48, 0xD5, 0xB3, 0x94, 0x9A, \
        0x18, 0x5B, 0x47, 0x59, 0xB1, 0x62, 0x04, 0x32, 0x33, 0x89, 0xE0, 0xCB, \
        0xF7, 0x33, 0xE1, 0x75, 0x0C, 0x22, 0xAC, 0x58, 0x93, 0x29, 0xE6, 0x6B, \
        0x00, 0xB0, 0x35, 0x07, 0x5B, 0xDD, 0xBE, 0x4C, 0x0A, 0x0B, 0xE7, 0x0B, \
        0xC9, 0xF0, 0x02, 0xB4, 0xB0, 0x7A, 0x42, 0x08 \
    };
#define BCH_MSG_MASK_M1                     0xF0
#define BCH_MSG_MASK_M2                     0x0F

// # Workspace memory requirements for the fast multiplication algorithms
//
// Karatsuba for multiplying two n-bit numbers (n > 1) requires memory
//   K(2) = 4
//   K(n) = 4 * ceil(n/2) + K(ceil(n/2))
// which solves to (for the worst case when n = 2^k+1):
//   K(n) <= 4*n - 4*k - 1.
// so that (note: TC_CUTOFF = 9)
//   K(9) <= 44,
// and actually equality is achieved here.
//
// Toom-Cook-3 (word aligned) for multiplying two n-bit numbers (n > TC2_CUTOFF) requires memory
//   T(10) = 54
//   T(n) = 11 * ceil(n/3) + 10 + T(ceil(n/3) + 2)
// which solves to (for the worst case when n = 2*3^k+4):
//   T(n) <= 11 * floor(n/2) + 32*k - 33
// so that (note: N_WORDS = 1108)
//   T(1108) <= 6253,
// although in reality T(1108) = 6177, which could be hardcoded since N_WORDS is fixed.
// We need space K(9) + T(1108):
#define GF2X_WORKSPACE_WORDS (44 + 6177)

#endif
