/**
 * test_bn.c - Unit tests for bn_add and bn_sub
 *
 * Compile with:
 *   gcc -I../include -o test_bn bn.c test_bn.c -Wall -Wextra
 * Run: ./test_bn
 */
#include <stdio.h>
#include <assert.h>
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

    printf("All tests passed!\n");
    return 0;
}
