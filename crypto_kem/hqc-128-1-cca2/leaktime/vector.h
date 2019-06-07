/**
 * \file vector.h
 * \brief Header file for vector.cpp
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 
#include <string.h>
#include <inttypes.h>

#include "parameters.h"
#include "rng.h"

#define BITMASK(a, size) ((1UL << (a % size)) - 1) /*!< Create a mask*/

void vect_fixed_weight(uint8_t* v, const uint16_t weight, AES_XOF_struct* ctx);
void vect_set_random(uint8_t* v, AES_XOF_struct* ctx);
void vect_set_random_from_randombytes(uint8_t* v);

void vect_add(uint8_t* o, uint8_t* v1, uint8_t* v2, uint32_t size);
int vect_compare(uint8_t* v1, uint8_t* v2, uint32_t size);

void vect_resize(uint8_t* o, uint32_t size_o, uint8_t* v, uint32_t size_v);
void vect_print(uint8_t* v, uint32_t size);

#endif