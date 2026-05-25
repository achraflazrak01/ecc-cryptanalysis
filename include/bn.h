/**
 * bn.h - Big Number arithmetic for 256-bit integers using 64-bit limbs
 * 
 * This is part of the ECC implementation (ECDH, ECDSA).
 * We use a fixed size of 4 limbs = 256 bits.
 * All functions are constant-time where needed.
 */

#ifndef BN_H
#define BN_H

#include <stdint.h> // for uint64_t, uint32_t, etc.
#include <stddef.h> //  for size_t
#include <stdbool.h> // for bool, true, false

/* ------------------------------------------------------------------------
 * Configuration
 * ------------------------------------------------------------------------ */

 // Number of 64-bit limbs to represent a 256-bit integer
 // 256 / 64 = 4
 #define BN_NLIMBS 4

// A big integer: array of limbs.
// Limbs are stored in little-endian order: limb[0] = least significant.
typedef uint64_t bn_limb_t;

// The main big integer type
typedef struct {
    bn_limb_t limbs[BN_NLIMBS];
} bn_t;

/* ------------------------------------------------------------------------
 * Constant-time helpers (inline for speed and clarity)
 * ------------------------------------------------------------------------ */

 // Return a mask that is all 1's if 'x' == 0, else all 0's.
 // This uses a trick: (x == 0) -> 1, then subtract 1 -> 0, than negate -> all ones.
 // For x != 0 -> (1) -> subtract 1 -> 0 -> negate -> 0.
 static inline uint64_t bn_is_zero_mask(uint64_t x) {
    return - (x == 0);
 }

 // Conditional move; if 'cond' is 1, copy 'src' into 'dst'; else keep 'dst' unchanged.
 // 'cond' must be 0 or 1. Works on whole bn_t struct via XOR trick.
 static inline void bn_cmov(bn_t *dst, const bn_t *src, uint64_t cond) {
    // 'cond' = 0 or 1. Create mask = 0 if cond=0, all ones if cond=1.
    uint64_t mask = -cond; // two's complement: -1 == all bits 1
    for (int i = 0; i < BN_NLIMBS; i++)
    {
        dst->limbs[i] ^= mask & (dst->limbs[i] ^ src->limbs[i]);
    }
 }

/* ------------------------------------------------------------------------
 * Basic arithmetic prototypes (implemented in bn.c)
 * ------------------------------------------------------------------------ */

// Set a big number to zero
void bn_zero(bn_t *a);

// Set from a native 64-bit integer (small value)
void bn_from_u64(bn_t *a, uint64_t x);

// Copy b into a
void bn_copy(bn_t *a, const bn_t *b);

// Compare two big numbers (constant-time) -> returns 1 if a > b, 0 if a == b, -1 if a < b
int bn_cmp(const bn_t *a, const bn_t *b);

// Addition a = b + c Return carry (0 or 1)
uint64_t bn_add(bn_t *a, const bn_t *b, const bn_t *c);

// Subtraction: a = b - c. Returns borrow (0 or 1) -> 0 = success, 1 = underflow
uint64_t bn_sub(bn_t *a, const bn_t *b, const bn_t *c);

// Greatest common divisor (Euclid) - for big integers
void bn_gcd(bn_t *g, const bn_t *a, const bn_t *b);

// Modular inversion using extended Euclidean algorithm
void bn_inv_mod(bn_t *inv, const bn_t *a, const bn_t *mod);

/* ------------------------------------------------------------------------
 * Helpers for binary GCD
 * ------------------------------------------------------------------------ */

 // Test if a big number is zero
 bool bn_is_zero(const bn_t *a);

 // Test if a big number is even (lowest bit == 0)
 bool bn_is_even(const bn_t *a);

 // In-place right shift by one bit (divide by 2)
 void bn_rshift1(bn_t *a);

 // In-place left shift by one bit (multiply by 2)
 void bn_rshift1(bn_t *a);

/* ------------------------------------------------------------------------
 * GCD (binary Euclidean algorithm)
 * ------------------------------------------------------------------------ */
// Compute g = gcd(a, b)
void bn_gcd(bn_t *g, const bn_t *a, const bn_t *b);

#endif /* BH_H */