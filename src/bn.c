/**
 * bn.c - Big number arithmetic for fixed-size 256-bit integers (4 limbs)
 *
 * Implements addition, subtraction, and basic utilities.
 * Carry/borrow propagation follows schoolbook algorithm.
 */

 #include "bn.h"
 #include <string.h> // for memset

/* ------------------------------------------------------------------------
 * Basic utilities
 * ------------------------------------------------------------------------ */

 // Set all bimbs to zero
 void bn_zero(bn_t *a) {
    for(int i = 0; i < BN_NLIMBS; i++) {
        a->limbs[i] = 0;
    }
 }

// Set from a smallnative integer (fits in 64 bits)
void bn_from_u64(bn_t *a, uint64_t x) {
    bn_zero(a);
    a->limbs[0] = x; // Only lowest limb set higher limbs remain 0
}

// Copy b into a
void bn_copy(bn_t *a, const bn_t *b) {
    for(int i = 0; i < BN_NLIMBS; i++) {
        a->limbs[i] = b->limbs[i];
    }
}

// Compare two big numbers.
// Returns: 1 if a > b, 0 if a == b,-1 if a < b
int bn_cmp(const bn_t *a, const bn_t *b) {
    // Compare from most signaficant limb down to least.
    for (int i = BN_NLIMBS - 1; i >= 0; i--)
    {
        if (a->limbs[i] > b->limbs[i]) return 1;
        if (a->limbs[i] < b->limbs[i]) return -1;
    }
    return 0;
    
}

/* ------------------------------------------------------------------------
 * Addition with carry
 * ------------------------------------------------------------------------
 * Computes a = b + c.
 * Returns carry out (0 or 1). The result a is stored in little-endian order.
 */
uint64_t bn_add(bn_t *a, const bn_t *b, const bn_t *c) {
    uint64_t carry = 0;
    for (uint64_t i = 0; i < BN_NLIMBS; i++) {
        // Use unsigned 128-bit intermediate to detect carry
        __uint128_t sum = (__uint128_t)b->limbs[i] + c->limbs[i] + carry;
        a->limbs[i] = (uint64_t)sum;          // lower 64 bits
        carry = (uint64_t)(sum >> 64);        // higher bits become new carry
    }
    return carry;  // 0 or 1 (since 2*2^64 fits in 65 bits)
    
}

/* ------------------------------------------------------------------------
 * Subtraction with borrow
 * ------------------------------------------------------------------------
 * Computes a = b - c.
 * Returns borrow out (0 if no borrow, 1 if underflow occurred).
 * If borrow = 1, the result is mathematically b - c + 2^(256).
 */
uint64_t bn_sub(bn_t *a, const bn_t *b, const bn_t *c) {
    uint64_t borrow = 0;
    for (int i = 0; i < BN_NLIMBS; i++) {
        // Use signed 128-bit to detect borrow: diff = b_i - c_i - borrow
        __int128_t diff = (__int128_t)b->limbs[i] - c->limbs[i] - borrow;
        if (diff < 0) {
            a->limbs[i] = (uint64_t) (diff + ((__int128_t)1 << 64));
            borrow = 1;
        } else {
            a->limbs[i] = (uint64_t)diff;
            borrow = 0;
        }
    }
    return borrow;
}