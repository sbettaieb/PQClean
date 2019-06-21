#ifndef BCH_H
#define BCH_H

/**
 * \file bch.h
 * \brief Header file for bch.c
 */

#include <stdint.h>

#include "parameters.h"

/**
 * \struct gf_tables
 * \brief a structure of a Galois Field Log and Anti-Log tables
 *
 * This structure allows to storage Log and Anti-Log tables.
 */
typedef struct gf_tables {
    uint16_t size; /*!< The size of the arrays of this structure*/
    int16_t log_tab[PARAM_GF_MUL_ORDER + 1]; /*!< An array that contains the Log values */
    int16_t antilog_tab[PARAM_GF_MUL_ORDER + 1]; /*!< An array that contains the Anti-Log values */
} gf_tables;

/**
 * \struct sigma_poly
 * \brief a structure of a Galois Field polynomial
 *
 * This structure allows to storage of a polynomial with coordinates in \f$ GF(2^{10}) \f$. We use
 * tis structure to compute the error location polynomial in the decoding phase of BCH code.
 */
typedef struct sigma_poly {
    uint16_t dim; /*!< The size of the array value of this structure*/
    uint16_t deg; /*!< The degree of the polynomial stored in the array value*/
    int16_t value[2 * PARAM_DELTA]; /*!< An array that contains the coordinates of the polynomial*/
} sigma_poly;

/**
 * \struct syndrome_set
 * \brief a structure of a Syndromes set
 *
 * This structure allows to storage of a set of syndromes.
 */
typedef struct syndrome_set {
    uint16_t size; /*!< The size of the array tab of this structure*/
    int16_t tab[2 * PARAM_DELTA]; /*!< An array that contains the values of syndromes*/
} syndrome_set;

/**
 * \struct cyclotomic_sets
 * \brief a structure of Cyclotomic sets
 *
 * This structure allows the storage of an array of cyclotomic sets
 */
typedef struct cyclotomic_sets {
    int8_t nb_c ;
    int8_t tab[PARAM_GF_MUL_ORDER + 1];
} cyclotomic_sets;

int16_t gf_get_antilog(gf_tables *tables, int16_t i);

int16_t gf_get_log(gf_tables *tables, int16_t i);

int16_t gf_mod(int16_t i);

void bch_code_encode(uint8_t *em, uint8_t *m);

void message_to_array(uint8_t *o, uint8_t *v);

void lfsr_encoder(uint8_t *em, uint8_t *g, uint8_t *m);

void array_to_codeword(uint8_t *v, uint8_t *c);

void bch_code_decode(uint8_t *m, uint8_t *em);

void gf_tables_init(gf_tables *tables);

void gf_generation(gf_tables *gf_tables);

void syndrome_init(syndrome_set *synd_set);

void syndrome_gen(syndrome_set *synd_set, gf_tables *tables, uint8_t *v);

void sigma_poly_init(sigma_poly *poly);

void get_error_location_poly(sigma_poly *sigma, gf_tables *tables, syndrome_set *synd_set);

void sigma_poly_copy(sigma_poly *p1, sigma_poly *p2);

void chien_search(uint16_t *error_pos, uint16_t *size, gf_tables *tables, sigma_poly *sigma);

void error_poly_gen(uint8_t *e, uint16_t *error_pos, uint16_t size);

void message_from_codeword(uint8_t *o, uint8_t *v);

void gf_tables_clear(gf_tables *gf_tables);

void syndrome_clear(syndrome_set *synd_set);

void sigma_poly_clear(sigma_poly *poly);

int16_t gf_mult(gf_tables *tables, int16_t a, int16_t b);

void cyclotomic_init(cyclotomic_sets *c_tab);

void cyclotomic_clear(cyclotomic_sets *c_tab);

void cyclotomic_gen(cyclotomic_sets *c_tab);

void compute_generator_poly(int16_t *g);

#endif
