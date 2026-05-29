/**
 * @file bn.h
 * @brief Big integer arithmetic for 256-bit fixed-size integers (4 limbs).
 * @details Implements addition, subtraction, comparison, GCD, modular inverse.
 *          All functions are constant-time where noted.
 * 
 * Limb representation: little-endian uint64_t[BN_NLIMBS]. 
 * BN_NLIMBS = 4 => 256 bits.
 * 
 * @author Achraf Lazrak
 * @date 2026
 */

#ifndef BN_H
#define BN_H

#include <stdint.h> ///> for uint64_t, uint32_t, etc.
#include <stddef.h> ///>  for size_t
#include <stdbool.h> ///> for bool, true, false

/* -------------------------- Configuration ------------------------- */

 #define BN_NLIMBS 4  ///< Number of 64-bit limbs to represent a 256-bit integer

typedef uint64_t bn_limb_t; ///< Limb type

/**
 * @brief Big integer type (fixed size).
 */
typedef struct {
    bn_limb_t limbs[BN_NLIMBS];  ///< Little-endian limbs
} bn_t;

/* ------ Constant-time helpers (inline for speed and clarity) ------ */

/**
 * @brief Constant-time test if a 64-bit word is zero.
 * @param x Input word.
 * @return Mask: all ones if x==0, else zero.
 */
 static inline uint64_t bn_is_zero_mask(uint64_t x) {
    return - (x == 0);
 }

/**
 * @brief Constant-time conditional move.
 * @param dst Destination.
 * @param src Source.
 * @param cond Must be 0 or 1. If 1, copy src to dst; else keep dst.
 */
 static inline void bn_cmov(bn_t *dst, const bn_t *src, uint64_t cond) {
    // 'cond' = 0 or 1. Create mask = 0 if cond=0, all ones if cond=1.
    uint64_t mask = -cond; ///> two's complement: -1 == all bits 1
    for (int i = 0; i < BN_NLIMBS; i++)
    {
        dst->limbs[i] ^= mask & (dst->limbs[i] ^ src->limbs[i]);
    }
 }

/* ------------------- Initialization / Utilities ------------------- */

/* Set a big number to zero. */
void bn_zero(bn_t *a);

/* Set from a native 64-bit integer (small value). */
void bn_from_u64(bn_t *a, uint64_t x);

/* Copy b into a. */
void bn_copy(bn_t *a, const bn_t *b);

/**
 * @brief Compare two big numbers.
 * @return 1 if a > b, 0 if equal, -1 if a < b.
 */
int bn_cmp(const bn_t *a, const bn_t *b);

 /* Test if a == 0. */
 bool bn_is_zero(const bn_t *a);

 /* Test if a big number is even (lowest bit = 0). */
 bool bn_is_even(const bn_t *a);

 /* Right shift by one bit (divide by 2). */
 void bn_rshift1(bn_t *a);

 /* Left shift by one bit (multiply by 2). */
 void bn_lshift1(bn_t *a);

/* ------------------- Basic Arithmetic ------------------- */

/**
 * @brief Addition: a = b + c.
 * @return Carry out (0 or 1).
 */
uint64_t bn_add(bn_t *a, const bn_t *b, const bn_t *c);

/**
 * @brief Subtraction: a = b - c.
 * @return Borrow out (0 if no borrow, 1 if underflow).
 */
uint64_t bn_sub(bn_t *a, const bn_t *b, const bn_t *c);

/* ------------------- GCD and Modular Inverse ------------------- */

/**
 * @brief Compute greatest common divisor: g = gcd(a, b).
 * @param g Output GCD.
 * @param a First operand.
 * @param b Second operand.
 * @note Uses binary (Stein) algorithm. Works for zero inputs.
 */
void bn_gcd(bn_t *g, const bn_t *a, const bn_t *b);

/**
 * @brief Compute modular inverse: inv = a^{-1} mod mod.
 * @param inv Output inverse (set only if return true).
 * @param a Input integer (must be reduced or not).
 * @param mod Modulus (must be odd and > 1).
 * @return true if inverse exists (gcd(a, mod) == 1), else false.
 * @note Uses binary extended Euclidean algorithm. 
 *       Runs in O(log^2 n). Not constant-time.
 */
bool bn_inv_mod(bn_t *inv, const bn_t *a, const bn_t *mod);

#endif /* BH_H */