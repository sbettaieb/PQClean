/**
 * \file gf2x.h
 * \brief Header file for gf2x.cpp
 */

#ifndef GF2X_H
#define GF2X_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#include <cinttypes>
#else
#define EXTERNC
#include <inttypes.h>
#endif

EXTERNC void ntl_cyclic_product(uint8_t*o, const uint8_t* v1, const uint8_t* v2);

#undef EXTERNC

#endif
