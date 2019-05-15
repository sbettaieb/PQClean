#ifndef SPX_THASH_H
#define SPX_THASH_H

#include <stdint.h>

void PQCLEAN_SPHINCSSHA256192FSIMPLE_CLEAN_thash_1(
    uint8_t *seed_state,
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]);

void PQCLEAN_SPHINCSSHA256192FSIMPLE_CLEAN_thash_2(
    uint8_t *seed_state,
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]);

void PQCLEAN_SPHINCSSHA256192FSIMPLE_CLEAN_thash_WOTS_LEN(
    uint8_t *seed_state,
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]);

void PQCLEAN_SPHINCSSHA256192FSIMPLE_CLEAN_thash_FORS_TREES(
    uint8_t *seed_state,
    unsigned char *out, const unsigned char *in,
    const unsigned char *pub_seed, uint32_t addr[8]);

#endif
