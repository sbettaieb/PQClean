/**
 * \file gf2x.c
 * \brief Implementation of multiplication of two polynomials
 *
 * The algorithms in this file are based on the following resources:
 *
 * - [bgtz]: Richard P. Brent, Pierrick Gaudry, Emmanuel Thomé, Paul Zimmermann.
 *           Faster Multiplication in GF(2)[x]. ANTS 2008.
 * - [bod]:  Marco Bodrato. Towards Optimal Toom-Cook Multiplication for Univariate
 *           and Multivariate Polynomials in Characteristic 2 and 0. WAIFI 2007.
 */
// FIXME There are implicit assumptions about the multiplicants being balanced.
// The methods appear to work even when a_len == b_len + 2 or a_len == b_len + 1,
// but this has not been rigourously tested.

#include "gf2x.h"

#include "parameters.h"

#include <assert.h>
#include <stdint.h>
#include <string.h>

typedef uint64_t word_t;

#define WORD_BITS (sizeof(word_t) << 3)
#define N_WORDS ((PARAM_N + WORD_BITS - 1) / WORD_BITS)
#define TAIL_BITS (PARAM_N % WORD_BITS)
#define TAIL_MASK ((((word_t)1) << TAIL_BITS) - 1)
#define PRE_BITS 4
#define PRE_MASK ((1 << PRE_BITS) - 1)
// REPAIR_MASK is hardcoded for (WORD_BITS, PRE_BITS) = (64, 4)
#define REPAIR_MASK 0x7777777777777777
#define TC2_CUTOFF 9

static void zero(word_t *c, size_t c_len);
static void assign(word_t *c, const word_t *a, size_t a_len);
static void add(word_t *c, const word_t *a, size_t a_len, const word_t *b, size_t b_len);
static void add_inplace(word_t *c, const word_t *a, size_t a_len);
static void mul_tc1(word_t *c, const word_t *a, size_t a_len, word_t b);
static void mul_tc2(word_t *c, const word_t *a, size_t a_len, const word_t *b, size_t b_len);
static void div_1_plus_W(word_t *c, size_t n);
static void mul_tc3w(word_t *c, const word_t *a, size_t a_len, const word_t *b, size_t b_len);
static void mul(word_t *c, const word_t *a, size_t a_len, const word_t *b, size_t b_len);
static void reduce(word_t *c);

// c = 0
static void zero(word_t *c, size_t c_len) {
    memset(c, 0, c_len * sizeof(word_t));
}

// c = a
static void assign(word_t *c, const word_t *a, size_t a_len) {
    memcpy(c, a, a_len * sizeof(word_t));
}

// c = a ^ b
static void add(word_t *c,
                const word_t *a, size_t a_len,
                const word_t *b, size_t b_len) {
    assert(c <= a || a + a_len <= c);
    assert(c <= b || b + b_len <= c);

    if (a_len < b_len) {
        add(c, b, b_len, a, a_len);
        return;
    }
    for (size_t i = 0; i < b_len; ++i) {
        c[i] = a[i] ^ b[i];
    }
    for (size_t i = b_len; i < a_len; ++i) {
        c[i] = a[i];
    }
}

// c ^= a
static void add_inplace(word_t *c, const word_t *a, size_t a_len) {
    assert(c <= a || a + a_len <= c);
    for (size_t i = 0; i < a_len; ++i) {
        c[i] ^= a[i];
    }
}

// c := a * b
// Extension of mul1 [bgtz], multiply n-by-1 words.
// Precomputes only once.
// c: (n+1)-word output
// a: n-word input
// b: one word input
static void mul_tc1(word_t *c, const word_t *a, size_t a_len, word_t b) {
    word_t table[1 << PRE_BITS] = {0, b};

    assert(c <= a + 1 || a + a_len <= c);

    memset(c, 0, (a_len + 1) * sizeof(word_t));

    // precompute
    for (size_t i = 2; i < (1 << PRE_BITS); i += 2) {
        table[i] = 2 * table[i / 2];
        table[i + 1] = table[i] ^ b;
    }

    for (size_t j = 0; j < a_len; ++j) {
        word_t aj = a[j];

        // multiply
        c[j] ^= table[aj & PRE_MASK];
        for (size_t i = PRE_BITS; i < WORD_BITS; i += PRE_BITS) {
            word_t tmp = table[(aj >> i) & PRE_MASK];
            c[j] ^= tmp << i;
            c[j + 1] ^= tmp >> (WORD_BITS - i);
        }

        // repair (note the bug in [bgtz])
        for (size_t i = 1; i < PRE_BITS; ++i) {
            aj = (aj >> 1) & REPAIR_MASK;
            c[j + 1] ^= aj & (-((b >> (WORD_BITS - i)) & 1));
        }
    }
}

// c := a * b
// Few-word multiplication using Karatsuba (Toom-2)
// TODO: assert no overlapping memory
static void mul_tc2(word_t *c,
                    const word_t *a, size_t a_len,
                    const word_t *b, size_t b_len) {
    // FIXME: memory overhead
    word_t a01[TC2_CUTOFF];
    word_t b01[TC2_CUTOFF];
    word_t cc[2 * TC2_CUTOFF];
    size_t n = (a_len + 1) / 2;

    add(a01, a, n, a + n, a_len - n);
    add(b01, b, n, b + n, b_len - n);
    mul(cc, a01, n, b01, n);
    mul(c, a, n, b, n);
    mul(c + 2 * n, a + n, a_len - n, b + n, b_len - n);
    add_inplace(cc, c, 2 * n);
    add_inplace(cc, c + 2 * n, a_len + b_len - 2 * n);
    add_inplace(c + n, cc, 2 * n);
}

// c := c / (1 + x**WORD_BITS)
static void div_1_plus_W(word_t *c, size_t n) {
    for (size_t i = 1; i < n - 1; ++i) {
        c[i] ^= c[i - 1];
    }
    // check for zero remainder
    assert(c[n - 2] == c[n - 1]);
}

// c := a * b
// Many-word multiplication using Toom-Cook-3
// Based on TC3W [bgtz] and balanced Toom-3 [bod]
static void mul_tc3w(word_t *c,
                     const word_t *a, size_t a_len,
                     const word_t *b, size_t b_len) {
    // FIXME: too much memory overhead (VLA is not an option)
    word_t c0[N_WORDS], c1[N_WORDS], c2[N_WORDS], c3[N_WORDS], c4[N_WORDS], c5[N_WORDS];
    size_t n = (a_len + 2) / 3;
    size_t c4_len;
    // a == a0 + a1X + a2X**2
    // b == b0 + b1X + b2X**2

    // notation: result [size] = ...

    // c0 [n+2] = a1*W + a2*W**2
    c0[0] = 0;
    assign(c0 + 1, a + n, n);
    c0[n + 1] = 0;
    add_inplace(c0 + 2, a + 2 * n, a_len - 2 * n);
    // c4 [n+2] = b1*W + b2*W**2
    c4[0] = 0;
    assign(c4 + 1, b + n, n);
    c4[n + 1] = 0;
    add_inplace(c4 + 2, b + 2 * n, b_len - 2 * n);
    // c5 [n] = a0 + a1 + a2
    add(c5, a, n, a + n, n);
    add_inplace(c5, a + 2 * n, a_len - 2 * n);
    // c2 [n] = b0 + b1 + b2
    add(c2, b, n, b + n, n);
    add_inplace(c2, b + 2 * n, b_len - 2 * n);
    // c1 [2n] = c2 * c5
    mul(c1, c2, n, c5, n);
    // c5 [n+2] = c5 + c0
    c5[n] = 0;
    c5[n + 1] = 0;
    add_inplace(c5, c0, n + 2);
    // c2 [n+2] = c2 + c4
    c2[n] = 0;
    c2[n + 1] = 0;
    add_inplace(c2, c4, n + 2);
    // c0 [n+2] = c0 + a0
    add_inplace(c0, a, n);
    // c4 [n+2] = c4 + b0
    add_inplace(c4, b, n);
    // c3 [2n+4] = c2 * c5
    mul(c3, c2, n + 2, c5, n + 2);
    // c2 [2n+4] = c0 * c4
    mul(c2, c0, n + 2, c4, n + 2);
    // c0 [2n] = a0 * b0
    mul(c0, a, n, b, n);
    // c4 [c4_len] = a2 * b2
    mul(c4, a + 2 * n, a_len - 2 * n, b + 2 * n, b_len - 2 * n);
    c4_len = a_len + b_len - 4 * n;
    // c3 [2n+2] = c3 + c2
    assert(c3[2 * n + 2] == c2[2 * n + 2] && c3[2 * n + 3] == c2[2 * n + 3]);
    add_inplace(c3, c2, 2 * n + 2);
    // c2 [2n+4] = c2 + c0
    add_inplace(c2, c0, 2 * n);
    // c2 [2n+3] = c2 / W + c3
    assert(c2[0] == 0); // c2 is now at c2+1
    add_inplace(c2 + 1, c3, 2 * n + 2);
    // c2 [2n] = (c2 + (1 + W**3)c4) / (1 + W)
    assert((c4_len < 2 * n - 1 && (c2 + 1)[2 * n + 1] == 0) || (c2 + 1)[2 * n + 1] == c4[2 * n - 2]);
    assert((c4_len < 2 * n   && (c2 + 1)[2 * n + 2] == 0) || (c2 + 1)[2 * n + 2] == c4[2 * n - 1]);
    add_inplace(c2 + 4, c4, c4_len);
    add_inplace(c2 + 1, c4, c4_len);
    div_1_plus_W(c2 + 1, 2 * n + 1);
    // c1 [2n] = c1 + c0
    add_inplace(c1, c0, 2 * n);
    // c3 [2n+2] = c3 + c1
    add_inplace(c3, c1, 2 * n);
    // c3 [2n] = c3 / (W**2 + W)
    assert(c3[0] == 0); // c3 is now at c3+1
    div_1_plus_W(c3 + 1, 2 * n + 1);
    // c1 [2n] = c1 + c2 + c4
    add_inplace(c1, c2 + 1, 2 * n);
    add_inplace(c1, c4, c4_len);
    // c2 [2n] = c2 + c3
    add_inplace(c2 + 1, c3 + 1, 2 * n);
    // c = c0 + c1*X + c2*X**2 + c3*X**3 + c4*X**4
    assign(c, c0, 2 * n);
    assign(c + 2 * n, c2 + 1, 2 * n);
    assign(c + 4 * n, c4, c4_len);
    add_inplace(c + n, c1, 2 * n);
    add_inplace(c + 3 * n, c3 + 1, 2 * n);
}

// generic multiplication function: depending on the sizes of the
// input it will call the fastest function for that size.
static void mul(word_t *c, const word_t *a, size_t a_len, const word_t *b, size_t b_len) {
    if (b_len == 0) {
        zero(c, a_len); // TODO: inefficient, better have a separate subroutine so this doesn't happen.
    } else if (b_len == 1) {
        mul_tc1(c, a, a_len, b[0]);
    } else if (b_len < TC2_CUTOFF) {
        mul_tc2(c, a, a_len, b, b_len);
    } else {
        mul_tc3w(c, a, a_len, b, b_len);
    }
}

// c := c mod (x**PARAM_N - 1)
// c can have up to twice as many words as required for x**(PARAM_N-1)
static void reduce(word_t *c) {
    c[0] ^= c[N_WORDS - 1] >> TAIL_BITS;
    c[N_WORDS - 1] &= TAIL_MASK;
    c[0] ^= c[N_WORDS] << (WORD_BITS - TAIL_BITS);
    c[1] ^= c[N_WORDS] >> TAIL_BITS;
    c[N_WORDS] = 0;
    for (size_t i = N_WORDS; i < 2 * N_WORDS; ++i) {
        c[i - N_WORDS] ^= c[i] << (WORD_BITS - TAIL_BITS);
        c[i + 1 - N_WORDS] ^= c[i] >> TAIL_BITS;
    }
    c[0] ^= c[N_WORDS - 1] >> TAIL_BITS;
    c[N_WORDS - 1] &= TAIL_MASK;
    c[0] ^= c[N_WORDS] << (WORD_BITS - TAIL_BITS);
}

/**
 * \fn void PQCLEAN_HQC1281CCA2_LEAKTIME_ntl_cyclic_product(uint8_t*o, const uint8_t* v1, const uint8_t* v2)
 * \brief Multiply two vectors
 *
 * Vector multiplication is defined as polynomial multiplication performed modulo the polynomial \f$ X^n - 1\f$.
 *
 * \param[out] o Product of <b>v1</b> and <b>v2</b>
 * \param[in] v1 Pointer to the first vector
 * \param[in] v2 Pointer to the second vector
 */
void PQCLEAN_HQC1281CCA2_LEAKTIME_ntl_cyclic_product(uint8_t *o, const uint8_t *v1, const uint8_t *v2) {
    word_t c[2 * N_WORDS], a[N_WORDS], b[N_WORDS];

    memcpy(a, v1, VEC_N_SIZE_BYTES);
    memcpy(b, v2, VEC_N_SIZE_BYTES);
    a[N_WORDS - 1] &= TAIL_MASK;
    b[N_WORDS - 1] &= TAIL_MASK;

    mul(c, a, N_WORDS, b, N_WORDS);
    reduce(c);

    memcpy(o, c, VEC_N_SIZE_BYTES);
}

