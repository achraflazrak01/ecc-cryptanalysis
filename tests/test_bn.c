/**
 * test_bn.c - Unit tests for bn_add and bn_sub
 *
 * Compile with:
 *   gcc -I../include -o test_bn bn.c test_bn.c -Wall -Wextra
 * Run: ./test_bn
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>   // for srand(), rand()
#include "bn.h"

// Helper to print a bn_t (for debugging)
void bn_print(const char *label, const bn_t *a) {
    printf("%s = ", label);
    for (int i = BN_NLIMBS - 1; i >= 0; i--) {
        printf("%016llX", a->limbs[i]);
        if (i > 0) printf(":");
    }
    printf("\n");
}

// Test two numbers: check that (b + c) produces expected value and carry
void test_add(const bn_t *b, const bn_t *c, const bn_t *expected, uint64_t expected_carry) {
    bn_t result;
    uint64_t carry = bn_add(&result, b, c);
    assert(carry == expected_carry);
    assert(bn_cmp(&result, expected) == 0);
}

// Test subtraction: b - c should equal expected, with borrow = expected_borrow
void test_sub(const bn_t *b, const bn_t *c, const bn_t *expected, uint64_t expected_borrow) {
    bn_t result;
    uint64_t borrow = bn_sub(&result, b, c);
    assert(borrow == expected_borrow);
    assert(bn_cmp(&result, expected) == 0);
}

// Simple 64-bit GCD for reference (naive Euclidean)
static uint64_t gcd64(uint64_t x, uint64_t y) {
    while (y != 0) {
        uint64_t t = y;
        y = x % y;
        x = t;
    }
    return x;
}

// Test bn_gcd against known cases and randoms pairs
void test_gcd(void) {
    bn_t a, b, g, expected;
    printf("Testing bn_gcd");

    // Test 1: gcd(0, 0) = 0
    bn_zero(&a);
    bn_zero(&b);
    bn_gcd(&g, &a, &b);
    assert(bn_is_zero(&g));

    // Test 2: gcd(0, 123) = 123
    bn_zero(&a);
    bn_from_u64(&b, 123);
    bn_gcd(&g, &a, &b);
    assert(bn_cmp(&g, &b) == 0);

    // Test 3: gcd(123, 0) = 123
    bn_gcd(&g, &b, &a);
    assert(bn_cmp(&g, &b) == 0);

    // Test 4: gcd(1, any) = 1
    bn_from_u64(&a, 1);
    bn_from_u64(&b, 0xFFFFFFFFFFFFFFFFULL);
    bn_gcd(&g, &a, &b);
    bn_from_u64(&expected, 1);
    assert(bn_cmp(&g, &expected) == 0);

    // Test 5: gcd(2^64, 2^64) = 2^64
    bn_zero(&a); a.limbs[1] = 1;   // 2^64
    bn_zero(&b); b.limbs[1] = 1;
    bn_gcd(&g, &a, &b);
    assert(bn_cmp(&g, &a) == 0);

    // Test 6: gcd(2^256-1, 0) = 2^256-1
    for (int i = 0; i < BN_NLIMBS; i++) a.limbs[i] = 0xFFFFFFFFFFFFFFFFULL;
    bn_zero(&b);
    bn_gcd(&g, &a, &b);
    assert(bn_cmp(&g, &a) == 0);

    // Test 7: gcd(2^256-1, 1) = 1
    bn_from_u64(&b, 1);
    bn_gcd(&g, &a, &b);
    bn_from_u64(&expected, 1);
    assert(bn_cmp(&g, &expected) == 0);

    // Test 8: gcd(6, 10) = 2 (low numbers)
    bn_from_u64(&a, 6);
    bn_from_u64(&b, 10);
    bn_gcd(&g, &a, &b);
    bn_from_u64(&expected, 2);
    assert(bn_cmp(&g, &expected) == 0);

    // Test 9: Random pairs within 64 bits (compare with gcd64)
    srand(123456); // fixed seed for reproducibility
    for (int test = 0; test < 100; test++) {
        uint64_t x = ((uint64_t)rand() << 32) | rand();
        uint64_t y = ((uint64_t)rand() << 32) | rand();
        if (x == 0 && y == 0) continue; // skip both zero
        
        bn_from_u64(&a, x);
        bn_from_u64(&b, y);
        bn_gcd(&g, &a, &b);

        uint64_t expected_gcd = gcd64(x, y);
        bn_from_u64(&expected, expected_gcd);
        assert(bn_cmp(&g, &expected) == 0);
    }

    // Test 10: Random 256-bit pairs - self-consistency:
    //   gcd(gcd(a,b), c) should equal gcd(a, gcd(b,c))
    bn_t c, g1, g2;
    for (int test = 0; test < 20; test++) {
        // Fill a,b,c with random limbs (use rand() for each limb)
        for (int i = 0; i < BN_NLIMBS; i++) {
            a.limbs[i] = ((uint64_t)rand() << 32) | rand();
            b.limbs[i] = ((uint64_t)rand() << 32) | rand();
            c.limbs[i] = ((uint64_t)rand() << 32) | rand();
        }
        // Avoid all zeros to keep gcd meaningful
        if (bn_is_zero(&a)) a.limbs[0] = 1;
        if (bn_is_zero(&b)) b.limbs[0] = 1;
        if (bn_is_zero(&c)) c.limbs[0] = 1;
        
        bn_gcd(&g1, &a, &b);
        bn_gcd(&g1, &g1, &c);   // g1 = gcd(gcd(a,b), c)

        bn_gcd(&g2, &b, &c);
        bn_gcd(&g2, &a, &g2);   // g2 = gcd(a, gcd(b,c))

        assert(bn_cmp(&g1, &g2) == 0);
    }
    printf("bn_gcd tests passed.\n");
}

// Test modular inverse: a * inv ≡ 1 (mod m)
void test_inv_mod(void) {
    bn_t a, mod, inv, one, expected;
    printf("Testing bn_inv_mod...\n");

    bn_from_u64(&one, 1);

    // Test 1: inv(3) mod 7 = 5, because 3*5=15=1 mod7
    bn_from_u64(&a, 3);
    bn_from_u64(&mod, 7);
    bool ok = bn_inv_mod(&inv, &a, &mod);
    assert(ok == true);
    bn_from_u64(&expected, 5);
    assert(bn_cmp(&inv, &expected) == 0);

    // Test 2: inv(1)
    bn_from_u64(&a, 1);
    bn_from_u64(&mod, 0xFFFFFFFFFFFFFFFFULL);
    ok = bn_inv_mod(&inv, &a, &mod);
    assert(ok == true);
    assert(bn_cmp(&inv, &a) == 0);

    // Test 3: non-invertiblz (gcd != 1) e.g. inv(2) mod 6
    bn_from_u64(&a, 2);
    bn_from_u64(&mod, 6);
    ok = bn_inv_mod(&inv, &a, &mod);
    assert(ok == false);

    // Test 4: random pairs (using small 64-bit modulus for which we can precompute inverse via extended Euclid)
    srand(12345);
    for (int test = 0; test < 100; test++) {
        uint64_t a_val = ((uint64_t)rand() << 32) | rand();
        uint64_t m_val = ((uint64_t)rand() << 32) | rand();
        if (m_val < 2) m_val = 2;
        if (m_val % 2) m_val++; // ensure odd for binary algorithm
        // Compute gcd64 and f gcd != 1 skip
        if (gcd64(a_val, m_val) != 1) continue;

        bn_from_u64(&a, a_val);
        bn_from_u64(&mod, m_val);
        ok = bn_inv_mod(&inv, &a, &mod);
        assert(ok == true);

        // Compute (a * inv) % m using 64-bit multiplication
        __uint128_t product = (__uint128_t)a_val * inv.limbs[0];
        uint64_t remainder = (uint64_t)(product % m_val);
        assert(remainder == 1);
    }

    bn_t p256;
    p256.limbs[0] = 0xFFFFFFFFFFFFFFFFULL;
    p256.limbs[1] = 0x00000000FFFFFFFFULL;  // because of the 2^224 term
    p256.limbs[2] = 0x0000000000000000ULL;
    p256.limbs[3] = 0xFFFFFFFF00000001ULL;

    bn_t two, p_plus_one, half;
    bn_from_u64(&two, 2);
    bn_copy(&p_plus_one, &p256);
    bn_from_u64(&one, 1);
    bn_add(&p_plus_one, &p_plus_one, &one);
    bn_copy(&half, &p_plus_one);
    bn_rshift1(&half);
    ok = bn_inv_mod(&inv, &two, &p256);
    assert(ok == true);
    assert(bn_cmp(&inv, &half) == 0);

    printf("bn_inv_mod tests passed.\n");

}

int main(void)
{
    bn_t b, c, expected;

    printf("Running bn_add / bn_sub unit tests...\n");

    // ---------- Addition tests ----------
    // Test 1: 0 + 0 = 0, carry 0
    bn_zero(&b);
    bn_zero(&c);
    bn_zero(&expected);
    test_add(&b, &c, &expected, 0);

    // Test 2: 1 + 1 = 2
    bn_from_u64(&b, 1);
    bn_from_u64(&c, 1);
    bn_from_u64(&expected, 2);
    test_add(&b, &c, &expected, 0);

    // Test 3: (2^64 - 1) + 1 = 2^64, i.e., lowest limb 0, limb[1] = 1, carry 0
    bn_from_u64(&b, 0xFFFFFFFFFFFFFFFFULL);
    bn_from_u64(&c, 1);
    bn_zero(&expected);
    expected.limbs[0] = 0;
    expected.limbs[1] = 1; // 2^64 fits into second limb
    test_add(&b, &c, &expected, 0);

    // Test 4: Ma value (all limbs 0xFFFF... ) + 1 = all limbs 0, carry 1
    for (int i = 0; i < BN_NLIMBS; i++) {
        b.limbs[i] = 0xFFFFFFFFFFFFFFFFULL;
        c.limbs[i] = 0;
    }
    c.limbs[0] = 1;
    bn_zero(&expected);   // result = 0, carry out = 1
    test_add(&b, &c, &expected, 1);

    // Test 5: Random 256-bit addition (simple: b = 0x123... , c = 0x456...)
    // We'll compute manually: 0x1111111111111111 + 0x2222222222222222 = 0x3333333333333333, no carry
    for (int i = 0; i < BN_NLIMBS; i++) {
        b.limbs[i] = 0x1111111111111111ULL;
        c.limbs[i] = 0x2222222222222222ULL;
        expected.limbs[i] = 0x3333333333333333ULL;
    }

    // ---------- Subtraction tests ----------
    // Test 6: 5 - 3 = 2, borrow 0
    bn_from_u64(&b, 5);
    bn_from_u64(&c, 3);
    bn_from_u64(&expected, 2);
    test_sub(&b, &c, &expected, 0);

    // Test 7: 3 - 5 = (3 - 5) mod 2^256 - 2, borrow 1
    bn_from_u64(&b, 3);
    bn_from_u64(&c, 5);
    // expected = 2^256 - 2, i.e., all lombs except lowest = 0xFFFFFFFFFFFFFFFE
    for (int i = 0; i < BN_NLIMBS; i++) {
        expected.limbs[i] = 0xFFFFFFFFFFFFFFFFULL;
    }
    expected.limbs[0] = 0xFFFFFFFFFFFFFFFEULL;
    test_sub(&b, &c, &expected, 1);

    // Test 8: 2^64 - 1 - 1 = 2^64 - 2, no borrow
    bn_from_u64(&b, 0xFFFFFFFFFFFFFFFFULL);
    bn_from_u64(&c, 1);
    bn_from_u64(&expected, 0xFFFFFFFFFFFFFFFEULL);
    test_sub(&b, &c, &expected, 0);

    // Test 9: 0 - 1 = all limbs oxFF...
    bn_zero(&b);
    bn_from_u64(&c, 1);
    for (int i = 0; i < BN_NLIMBS; i++) {
        expected.limbs[i] = 0xFFFFFFFFFFFFFFFFULL;
    }
    test_sub(&b, &c, &expected, 1);

    // Test 10: Equal numbers: 123 - 123 = 0, borrow 0
    bn_from_u64(&b, 123);
    bn_copy(&c, &b);
    bn_zero(&expected);
    test_sub(&b, &c, &expected, 0);

    test_gcd();

    printf("All tests passed!\n");
    return 0;
}
