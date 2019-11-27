/**
 * \file gf2x.c
 * \brief AVX2 implementation of multiplication of two polynomials
 */

#include "gf2x.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>

#define VECTORSIZE 256 /*!< The size in bits of a vector*/
#define BYTESPERVECTOR (VECTORSIZE/8) /*!< The number of bytes in a vector*/
#define BITSPERINT 32 /*!< The number of bits per int*/
#define INTSPERVECTOR (VECTORSIZE/BITSPERINT) /*!< The number of ints per vector*/
#define VEC_N_ARRAY_SIZE_VEC CEIL_DIVIDE(PARAM_N,256) /*!< The number of needed vectors to store PARAM_N bits*/
typedef uint32_t v8uint32 __attribute__ ((vector_size (BYTESPERVECTOR), aligned (BYTESPERVECTOR)));/*!< v8uint32: a vector of 8 unsigned ints*/
static void shiftXor(v8uint32 dest[], v8uint32 source[], uint32_t length, uint32_t distance);


/**
 * @brief Computes dest ^= source>>distance,
 *
 * Source and dest are bignumbers
 * stored little endian in memory and read as vectors
 * Routine is optimized to efficiently use vector registers
 * we assume source is copied twice, you can use it to rotate
 * in picture form (with four bit words):
 * say distance = 11 and, v2 has length 14 bits
 * v2:       aaaa bbbb cccc dd..
 * we made a double vector as make_double_vector
 * to do that, shiftxor will read three vectors and write:
 * <pre>
 * (0000,0000)aaaa,bbbb cccc,ddaa aabb,bbcc ccdd,aaaa bbbb,cccc
 *                         \ |||\ |||\ |||\ |
 * v2:                     aaaa bbbb cccc dd..
 * </pre>
 *
 * @param[out] dest Array that is the result
 * @param[in] source Array that is the source of data
 * @param[in] length Integer that is size of the source
 * @param[in] distance Integer that is the amount of shits to do
 */
static void shiftXor(v8uint32 dest[], v8uint32 source[], uint32_t length, uint32_t distance) {
    // Split the distance in number of bits, ints, and vectors
    // there are 32 bits in an int, and 8 ints in a vector
    const int32_t distBits = distance & 0x1f;
    const int32_t minDistBits = BITSPERINT - distBits;
    const int32_t distInts = ( distance >> 5 ) & 0x7;
    const int32_t distVectors = distance >> 8;
    // Prepare for shuffling integers in the vectors
    const v8uint32 countVector = {0, 1, 2, 3, 4, 5, 6, 7};
    const v8uint32 rightPerm = countVector + distInts;
    const v8uint32 leftPerm = countVector + distInts + 1;

    v8uint32 previousVector, currentVector, rightVector, leftVector, shifted;
    v8uint32 *sourcePtr = source + distVectors;
    v8uint32 *destPtr = dest;
    previousVector = *sourcePtr++;

    for (uint32_t offset = 0; offset < length; offset++) {
        currentVector = *sourcePtr++;
        // get the appropriate integers from the last two vectors
        rightVector = __builtin_shuffle(
                          previousVector, currentVector, rightPerm);
        leftVector = __builtin_shuffle(
                         previousVector, currentVector, leftPerm);
        // get the appropriate bits from the integers
        shifted = rightVector >> distBits | leftVector << minDistBits;
        // do the xor
        *destPtr++ ^= shifted;
        previousVector = currentVector;
    }
}



/**
 * @brief Construct vector concatenated with itself in aligned memory
 *
 * First we have to make two concatenated copies of the vector
 * in picture form, where we take:
 * PARAM_N = 13                 BYTESPERVECTOR = 2
 * VEC_N_ARRAY_SIZE_VEC = 2     VECTORSIZE = 8
 * VEC_N_ARRAY_SIZE_INT = 4
 * we have a number of 13 bits, written little-endian in words:
 * <pre>
 * v2:        aaaa bbbb cccc d...
 * we want to concatenate with itself to get:
 *            aaaa bbbb cccc daaa abbb bccc cd..
 * let's say a vector is two of these mini-ints, so we see this as:
 *            aaaa,bbbb cccc,d...
 * first we make room for twice the number of vectors, plus one:
 *  0000,0000 0000,0000 0000,0000 0000,0000 0000,0000
 * copy the data, skipping one initial vector:
 *  0000,0000 aaaa,bbbb cccc,d... 0000,0000 0000,0000
 *  now we use shiftXor with a distance of 2*8-13=3, and a length of 3:
 * source:   ^   dest: ^
 * 0000,0000 aaaa,bbbb cccc,d... 0000,0000 0000,0000
 *     ^ ^^^^ ^^^       vvvv vvvv      xor action 1
 *  0000,0000 aaaa,bbbb cccc,daaa 0000,0000 0000,0000
 *  xor action 2:^ ^^^^ ^^^       vvvv vvvv
 *  0000,0000 aaaa,bbbb cccc,daaa abbb,bccc 0000,0000
 *  xor action 3:          ^ ^^^^ ^^^       vvvv vvvv
 *  0000,0000 aaaa,bbbb cccc,daaa abbb,bccc cdaa,abbb
 *  (note the extra wraparound).
 *  The output is at result+1
 * </pre>
 *
 * @param[in] dv Pointer to the resulting double vector
 * @param[in] v Pointer to the input vector
 * \return a pointer to the double vector
 */
static inline void make_double_vector(v8uint32 *dv, const uint8_t *v) {
    memcpy(dv + 1, v, VEC_N_SIZE_BYTES);
    shiftXor(dv + VEC_N_ARRAY_SIZE_VEC, dv, VEC_N_ARRAY_SIZE_VEC + 1, VEC_N_ARRAY_SIZE_VEC * VECTORSIZE - PARAM_N);
}



/**
 * @brief Multiply two vectors
 *
 * This functions multiplies a sparse vector <b>v1</b> (of Hamming weight equal to <b>weight</b>)
 * and a dense vector <b>v2</b>. The multiplication is done using the matrix-vector form of the multiplication of two
 * polynomials modulo \f$ X^n - 1\f$. where the matrix is generated using the dense vector.
 *
 * @param[out] o Pointer to the result
 * @param[in] v1 Pointer to the first vector
 * @param[in] v2 Pointer to the first vector
 * @param[in] weight Integer that is the weigt of the sparse vector
 */
void PQCLEAN_HQC2561CCA2_AVX2_LEAKTIME_vect_mul(uint8_t *o, const uint32_t *v1, const uint8_t *v2, uint32_t weight) {
    v8uint32 double_v2 [2 * VEC_N_ARRAY_SIZE_VEC + 1] = {0};
    v8uint32 result [VEC_N_ARRAY_SIZE_VEC] = {0};

    make_double_vector(double_v2, v2);

    for (uint16_t i = 0 ; i < weight ; ++i) {
        shiftXor(result, double_v2 + 1, VEC_N_ARRAY_SIZE_VEC, PARAM_N - v1[i]);
    }

    memcpy(o, result, VEC_N_SIZE_BYTES);
    o[VEC_N_SIZE_BYTES - 1] &= BITMASK(PARAM_N, 8);
}
