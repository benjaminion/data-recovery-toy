/*
 * Copyright 2021 Benjamin Edgington
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This a toy example of data recovery as proposed for Ethereum 2.0, a la
 * https://ethresear.ch/t/reed-solomon-erasure-code-recovery-in-n-log-2-n-time-with-ffts/3039
 *
 * It is intended to be read alongside the accompanying writeup at https://hackmd.io/@benjaminion/data_recovery
 *
 * This does the data encoding over a finite field (integers mod 17), using 4 as the primitive fourth root of 1. This
 * gives the powers of roots of unity we use as [1, 4, 16, 13].
 *
 * Notes:
 *  - It is hardcoded to deal with only 4 samples (two data elements!), up to two of which can be lost.
 *  - The division is hardcoded for MOD = 17 because I am lazy.
 *
 * For a proper implementation, see my library https://github.com/benjaminion/c-kzg
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#define MOD 17

// Fourth root of 1
const int i = 4;

void print_poly(const char *s, const int *a) {
    printf("%18s: [%d, %d, %d, %d]\n", s, a[0], a[1], a[2], a[3]);
}

bool equal(int a, int b) {
    return a % MOD == b % MOD;
}

int neg(int a) {
    return (MOD - a);
}

int add(int a, int b) {
    return (a + b) % MOD;
}

int sub(int a, int b) {
    return (a + neg(b)) % MOD;
}

int mul(int a, int b) {
    return (a * b) % MOD;
}

int div(int a, int b) {
    int inv[] = {0, 1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16};
    int ret = (a * inv[b % MOD]) % MOD;
    assert(equal(a, mul(ret, b)));
    return ret;
}

// Forward Fourier Transform: converts polynomial coefficients into polynomial evaluations at roots of unity
// The roots of unity we are using are [1, i, -1, -i] in that order.
void eval_from_poly(int *eval, const int *coeffs) {
    int c0_p_c2 = add(coeffs[0], coeffs[2]);
    int c0_m_c2 = sub(coeffs[0], coeffs[2]);
    int c1_p_c3 = add(coeffs[1], coeffs[3]);
    int c1_m_c3 = sub(coeffs[1], coeffs[3]);
    eval[0] = add(c0_p_c2, c1_p_c3);
    eval[1] = add(c0_m_c2, mul(i, c1_m_c3));
    eval[2] = sub(c0_p_c2, c1_p_c3);
    eval[3] = sub(c0_m_c2, mul(i, c1_m_c3));
}

// Reverse Fourier Transform: converts polynomial evaluations at roots of unity into polynomial coefficients
// The roots of unity we are using are [1, i, -1, -i] in that order.
void poly_from_eval(int *coeffs, const int *eval) {
    int c0_p_c2 = add(eval[0], eval[2]);
    int c0_m_c2 = sub(eval[0], eval[2]);
    int c1_p_c3 = add(eval[1], eval[3]);
    int c1_m_c3 = sub(eval[1], eval[3]);
    coeffs[0] = add(c0_p_c2, c1_p_c3);
    coeffs[1] = sub(c0_m_c2, mul(i, c1_m_c3));
    coeffs[2] = sub(c0_p_c2, c1_p_c3);
    coeffs[3] = add(c0_m_c2, mul(i, c1_m_c3));
    for (int j = 0; j < 4; j++) coeffs[j] = div(coeffs[j], 4);
}

// Given polynomial p(x), return the polynomial p(x / k).
// Does this by multiplying the coefficients by powers of k.
void scale(int *a, int k) {
    int fac = k;
    a[1] = mul(a[1], fac);
    fac = mul(fac, k);
    a[2] = mul(a[2], fac);
    fac = mul(fac, k);
    a[3] = mul(a[3], fac);
}

// Given polynomial p(x), return the polynomial p(x * k).
// Does this by dividing the coefficients by powers of k.
void unscale(int *a, int k) {
    int fac = k;
    a[1] = div(a[1], fac);
    fac = mul(fac, k);
    a[2] = div(a[2], fac);
    fac = mul(fac, k);
    a[3] = div(a[3], fac);
}

int main() {
    // My data is [5, 7] - treat it as the the coefficients of a polynomial, D(x).
    // The data is two elements, and then extended with the same number of zeros.
    const int data_poly[] = {5, 7, 0, 0};
    print_poly("Initial values", data_poly);

    // Encode my data by evaluating the polynomial at the roots of unity.
    int data_eval[4];
    eval_from_poly(data_eval, data_poly);
    print_poly("Data encoded", data_eval);

    // Lose part of the data. In this case, lose
    // -  position 1, corresponding to r^1 = i
    // -  position 2, correspoding to r^2 = -1
    // But could be any two elements. r is our fourth-root of unity, i.
    data_eval[1] = 0;
    data_eval[2] = 0;
    print_poly("Data with missing", data_eval);

    // Construct the zero polynomial as the product of (1 - r^j) for each of the indices j that is missing.
    // Our zero poly is (x - r^1)(x - r^2) = (x - 4)(x - 16) = x^2 - 3x + 13 = x^2 + 14x + 13
    int zero_poly[] = {13, 14, 1, 0};
    int zero_poly_eval[4];
    eval_from_poly(zero_poly_eval, zero_poly);
    print_poly("ZeroPoly eval", zero_poly_eval);

    // Create the evaluation (E * Z)(r^j) for each j. It has zeros at indices where we lost data.
    int ez_eval[4];
    for (int j = 0; j < 4; j++) {
        ez_eval[j] = mul(data_eval[j], zero_poly_eval[j]);
    }
    print_poly("EZ eval", ez_eval);

    // Interpolate to get (E * Z)(x) = (D * Z)(x). The equality holds since, by construction, E and D agree both where
    // non-zero, and where zero due to Z.
    int dz_poly[4];
    poly_from_eval(dz_poly, ez_eval);
    print_poly("EZ = DZ poly", dz_poly);

    // Scale ("shift") the polynomials so that we can divide them without hitting a zero in the zero poly.
    // Any scale factor is ok, as long as it is not one of the roots of unity (or zero).
    scale(dz_poly, 2);
    scale(zero_poly, 2);
    print_poly("DZ poly scaled", dz_poly);
    print_poly("ZeroPoly scaled", zero_poly);

    // Now we will divide the scaled polynomial (D * Z)(x / 2) by Z(x / 2), which will result in D(x / 2) - our
    // (scaled) original data. We do this via convolution: convert to evaluation form, divide pointwise, and convert
    // back to polynomial form.

    // Convert (D * Z)(x / 2) and Z(x / 2) to evaluation form,
    int dz_scaled_eval[4];
    int zero_scaled_eval[4];
    eval_from_poly(dz_scaled_eval, dz_poly);
    eval_from_poly(zero_scaled_eval, zero_poly);
    print_poly("DZ eval scaled", dz_scaled_eval);
    print_poly("Zero eval scaled", zero_scaled_eval);

    // Divide pointwise
    int quotient_eval[4];
    for (int j = 0; j < 4; j++) {
        quotient_eval[j] = div(dz_scaled_eval[j], zero_scaled_eval[j]);
    }
    print_poly("Quotient eval", quotient_eval);

    // Convert back to polynomial form to get D(x / 2).
    int recovered_poly[4];
    poly_from_eval(recovered_poly, quotient_eval);
    print_poly("Scaled recovered", recovered_poly);

    // Reverse our earlier scaling to recover D(x).
    unscale(recovered_poly, 2);
    print_poly("Recovered values", recovered_poly);

    return 0;
}
