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

/* ------------------------------------------------------------------------
 * Binary GCD (Stein's algorithm)
 * ------------------------------------------------------------------------
 * Computes g = gcd(a, b). Uses only subtraction and halving.
 * The algorithm:
 *   - Remove common factors of 2 (d)
 *   - Then repeatedly: if both even -> half; if one even -> half; else subtract smaller from larger.
 *   - Finally multiply back by 2^d.
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

/* ------------------------------------------------------------------------
 * Binary extended GCD (modular inverse)
 * ------------------------------------------------------------------------
 * For odd modulus 'mod', compute inv such that (a * inv) % mod == 1.
 * Uses the binary (Stein) algorithm with co-factors.
 * Returns true on success (inverse exists), false otherwise.
 */
bool bn_inv_mod(bn_t *inv, const bn_t *a, const bn_t *mod) {
    bn_t u, v, x1, x2, t, q, tmp;
    int shift;

    // Copy arguments to work on
    bn_copy(&u, a);
    bn_copy(&v, mod);
    bn_from_u64(&x1, 1);
    bn_zero(&x2);

    // Reduce u modulo mod (ensure u < mod)
    while (bn_cmp(&u, &v) >= 0) {
        bn_sub(&u, &u, &v);
    }

    // Main loop: while u != 0
    while (!bn_is_zero(&u)) {
        // Remove factors of 2 from u
        shift = 0;
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
            shift++;
        }
        
        // Remove factors of 2 from v
        while (bn_is_even(&v)) {
            bn_rshift1(&v);
            if(bn_is_even(&x2)) {
                bn_rshift1(&x2);
            } else {
                bn_add(&tmp, &x2, mod);
                bn_rshift1(&tmp);
                bn_copy(&x2, &tmp);
            }
        }

        // Now both u and v are odd
        if (bn_cmp(&u, &v) >= 0) {
            bn_sub(&u, &u, &v);
            bn_sub(&x1, &x1, &x2);
        } else {
            bn_sub(&v, &v, &u);
            bn_sub(&x2, &x2, &x1);
        }
    }
    
    // At this point, v = gcd(a, mod). If v != 1, no inverse.
    if(!(v.limbs[0] == 1 && bn_is_zero(&v) == false &&
            (v.limbs[1] | v.limbs[2] | v.limbs[3]) == 0)) {
            return false;  // gcd != 1
    }
    // Inverse is x2 mod mod. Ensure positive.
    bn_copy(inv, &x2);
    // If inv negative? Our algorithm keeps x2 in [0, mod-1] because we add mod when needed.
    // But final x2 might be > mod, so reduce once.
    while (bn_cmp(inv, mod) >= 0) {
        bn_sub(inv, inv, mod);
    }
    while (inv->limbs[0] & 0x8000000000000000ULL) { // handle sign if we used signed
        bn_add(inv, inv, mod);
    }
    return true;
}