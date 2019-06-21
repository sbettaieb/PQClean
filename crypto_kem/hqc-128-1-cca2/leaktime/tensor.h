/**
 * \file tensor.h
 * \brief Header file for tensor.c
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <stdint.h>

void tensor_code_encode(uint8_t *em, uint8_t *m);
void tensor_code_decode(uint8_t *m, uint8_t *em);

#endif
