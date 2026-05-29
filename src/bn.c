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

/**
 * @brief Addition with carry propagation.
 * @param a Output: b + c.
 * @param b First addend.
 * @param c Second addend.
 * @return Carry out (0 or 1).
 * @note Uses 128-bit intermediate for portable carry detection.
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

/* ------------------------------------------------------------------------
 * GCD helpers: zero, even, shifts
 * ------------------------------------------------------------------------ */

 bool bn_is_zero(const bn_t *a) {
    for (int i = 0; i < BN_NLIMBS; i++) {
        if (a->limbs[i] != 0) return false;
    }
    return true;
 }

 bool bn_is_even(const bn_t *a) {
    return (a->limbs[0] & 1) == 0;
 }

 // a = a >> 1 (right shift one bit, no carry)
 void bn_rshift1(bn_t *a) {
    uint64_t carry = 0;
    for (int i = BN_NLIMBS - 1; i >= 0; i--) {
        uint64_t limb = a->limbs[i];
        a->limbs[i] = (limb >> 1) | (carry << 63);
        carry = limb & 1;
    }
 }

 // a = a << 1 (left shift one bit)
 void bn_lshift1(bn_t *a) {
    uint64_t carry = 0;
    for (int i = 0; i < BN_NLIMBS; i++)
    {
        uint64_t limb = a->limbs[i];
        a->limbs[i] = (limb << 1) | carry;
        carry = limb >> 63;
    }
    // Note: we ignore overflow beyond 256 bits (no need for GCD)
 }

/**
 * @brief Binary GCD (Stein's algorithm).
 * @details Removes common factors of 2, then uses subtraction and halving.
 *          Faster than Euclidean algorithm for large numbers.
 */
void bn_gcd(bn_t *g, const bn_t *a, const bn_t *b) {
    bn_t u, v;
    int shift = 0;

    bn_copy(&u, a);
    bn_copy(&v, b);

    // GCD(0, v) = v
    if (bn_is_zero(&u)) {
        bn_copy(g, &v);
        return;
    }
    if (bn_is_zero(&v)) {
        bn_copy(g, &u);
        return;
    }

    // Remove common factors of 2
    while (bn_is_even(&u) && bn_is_even(&v)) {
        bn_rshift1(&u);
        bn_rshift1(&v);
        shift++;
    }

    while (!bn_is_zero(&u)) {
        // Remove factors of 2 from u
        while(bn_is_even(&u)) {
            bn_rshift1(&u);
        }
        // Remove factors of 2 from v
        while(bn_is_even(&v)) {
            bn_rshift1(&v);
        }
        
        // Now both u and v are odd
        if (bn_cmp(&u, &v) >= 0) {   // u >= v
            bn_sub(&u, &u, &v);
        } else {
            bn_sub(&v, &v, &u);
        }
    }

    // Result is v * 2^shift
    bn_copy(g, &v);
    while (shift--) {
        bn_lshift1(g);
    }
}

/**
 * @brief Modular inverse using binary extended GCD.
 * @details Implements algorithm 4.5.2 from "Handbook of Applied Cryptography".
 *          Reduces temporary copies by reusing inv as x2 and limiting temporaries.
 * @param inv Output: a^{-1} mod mod.
 * @param a Input integer.
 * @param mod Odd modulus.
 * @return true if inverse exists.
 */
bool bn_inv_mod(bn_t *inv, const bn_t *a, const bn_t *mod) {
    bn_t u, v, x1, tmp;

    // Use u = a %  mod, v = mod
    bn_copy(&u, a);
    bn_copy(&v, mod);
    bn_from_u64(&x1, 1);
    bn_zero(inv);

    // Reduce u modulo mod (ensure u < mod)
    while (bn_cmp(&u, &v) >= 0) {
        bn_sub(&u, &u, &v);
    }

    // Main loop: while u != 0
    while (!bn_is_zero(&u)) {
        // Remove factors of 2 from u
        while (bn_is_even(&u)) {
            bn_rshift1(&u);
            // Update x1: if x1 is even, half it; else (x1 + mod)/2
            if (bn_is_even(&x1)) {
                bn_rshift1(&x1);
            } else {
                bn_add(&tmp, &x1, mod);
                bn_rshift1(&tmp);
                bn_copy(&x1, &tmp);
            }
        }
        
        // Remove factors of 2 from v
        while (bn_is_even(&v)) {
            bn_rshift1(&v);
            if(bn_is_even(inv)) {
                bn_rshift1(inv);
            } else {
                bn_add(&tmp, inv, mod);
                bn_rshift1(&tmp);
                bn_copy(inv, &tmp);
            }
        }

        // Now both u and v are odd
        if (bn_cmp(&u, &v) >= 0) {
            bn_sub(&u, &u, &v);
            bn_sub(&x1, &x1, inv);
        } else {
            bn_sub(&v, &v, &u);
            bn_sub(inv, inv, &x1);
        }
    }
    
    // At this point, v = gcd(a, mod). If v != 1, no inverse.
    if(!(v.limbs[0] == 1 && bn_is_zero(&v) == false &&
            (v.limbs[1] | v.limbs[2] | v.limbs[3]) == 0)) {
            return false;  // gcd != 1
    }

    // If inv negative? Our algorithm keeps x2 in [0, mod-1] because we add mod when needed.
    while (bn_cmp(inv, mod) >= 0) {
        bn_sub(inv, inv, mod);
    }
    // If negative (unlikely with our algorithm), add mod once
    while (inv->limbs[0] & 0x8000000000000000ULL) {
        bn_add(inv, inv, mod);
    }
    return true;
}