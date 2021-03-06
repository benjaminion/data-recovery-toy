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
 * Notes:
 *  - It is hardcoded to deal with only 4 samples (two data elements!), up to two of which can be lost.
 *  - This is not implemented over a finite field as per Eth2, but over the complex integers.
 *
 * For a proper implementation, see my library https://github.com/benjaminion/c-kzg
 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

typedef int number;

typedef struct {
    number re;
    number im;
} complex;

const complex zero = {0, 0};
const complex i = {0, 1};
const complex two = {2, 0};
const complex four = {4, 0};

void print_c(complex a) {
    printf("%d", a.re);
    if (a.im != 0) {
        if (a.im > 0) {
            if (a.im == 1)
                printf(" + i");
            else
                printf(" + %di", a.im);
        } else {
            if (a.im == -1)
                printf(" - i");
            else
                printf(" - %di", -a.im);
        }
    }
}

void print_poly(const char *s, const complex *a) {
    printf("%18s: [", s);
    print_c(a[0]);
    printf(", ");
    print_c(a[1]);
    printf(", ");
    print_c(a[2]);
    printf(", ");
    print_c(a[3]);
    printf("]\n");
}

bool equal(complex a, complex b) {
    return a.re == b.re && a.im == b.im;
}

complex new_complex(number re, number im) {
    complex ret = {re, im};
    return ret;
}

complex cconj(complex a) {
    return new_complex(a.re, -a.im);
}

complex add(complex a, complex b) {
    return new_complex(a.re + b.re, a.im + b.im);
}

complex sub(complex a, complex b) {
    return new_complex(a.re - b.re, a.im - b.im);
}

complex mul(complex a, complex b) {
    return new_complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

complex div(complex a, complex b) {
    complex b_conj = cconj(b);
    complex ab_conj = mul(a, b_conj);
    complex bb_conj = mul(b, b_conj);
    complex ret = new_complex(ab_conj.re / bb_conj.re, ab_conj.im / bb_conj.re);
    // Can fail due to the integer division, or division by zero:
    assert(equal(a, mul(ret, b)));
    return ret;
}

// Forward Fourier Transform: converts polynomial coefficients into polynomial evaluations at roots of unity
// The roots of unity we are using are [1, i, -1, -i] in that order.
void eval_from_poly(complex *eval, const complex *coeffs) {
    complex c0_p_c2 = add(coeffs[0], coeffs[2]);
    complex c0_m_c2 = sub(coeffs[0], coeffs[2]);
    complex c1_p_c3 = add(coeffs[1], coeffs[3]);
    complex c1_m_c3 = sub(coeffs[1], coeffs[3]);
    eval[0] = add(c0_p_c2, c1_p_c3);
    eval[1] = add(c0_m_c2, mul(i, c1_m_c3));
    eval[2] = sub(c0_p_c2, c1_p_c3);
    eval[3] = sub(c0_m_c2, mul(i, c1_m_c3));
}

// Reverse Fourier Transform: converts polynomial evaluations at roots of unity into polynomial coefficients
// The roots of unity we are using are [1, i, -1, -i] in that order.
void poly_from_eval(complex *coeffs, const complex *eval) {
    complex c0_p_c2 = add(eval[0], eval[2]);
    complex c0_m_c2 = sub(eval[0], eval[2]);
    complex c1_p_c3 = add(eval[1], eval[3]);
    complex c1_m_c3 = sub(eval[1], eval[3]);
    coeffs[0] = add(c0_p_c2, c1_p_c3);
    coeffs[1] = sub(c0_m_c2, mul(i, c1_m_c3));
    coeffs[2] = sub(c0_p_c2, c1_p_c3);
    coeffs[3] = add(c0_m_c2, mul(i, c1_m_c3));
    for (int j = 0; j < 4; j++) coeffs[j] = div(coeffs[j], four);
}

// Given polynomial p(x), return the polynomial p(x / k).
// Does this by multiplying the coefficients by powers of k.
void scale(complex *a, complex k) {
    complex fac = k;
    a[1] = mul(a[1], fac);
    fac = mul(fac, k);
    a[2] = mul(a[2], fac);
    fac = mul(fac, k);
    a[3] = mul(a[3], fac);
}

// Given polynomial p(x), return the polynomial p(x * k).
// Does this by dividing the coefficients by powers of k.
void unscale(complex *a, complex k) {
    complex fac = k;
    a[1] = div(a[1], fac);
    fac = mul(fac, k);
    a[2] = div(a[2], fac);
    fac = mul(fac, k);
    a[3] = div(a[3], fac);
}

int main() {
    // My data is [5, 7] - treat it as the the coefficients of a polynomial, D(x).
    // The data is two elements, and then extended with the same number of zeros.
    const complex data_poly[] = {{5, 0}, {7, 0}, {0, 0}, {0, 0}};
    print_poly("Initial values", data_poly);

    // Encode my data by evaluating the polynomial at the roots of unity.
    complex data_eval[4];
    eval_from_poly(data_eval, data_poly);
    print_poly("Data encoded", data_eval);

    // Lose part of the data. In this case, lose
    // -  position 1, corresponding to r^1 = i
    // -  position 2, correspoding to r^2 = -1
    // But could be any two elements. r is our fourth-root of unity, i.
    data_eval[1] = new_complex(0, 0);
    data_eval[2] = new_complex(0, 0);
    print_poly("Data with missing", data_eval);

    // Construct the zero polynomial as the product of (1 - r^j) for each of the indices j that is missing.
    // Our zero poly is (x - r^1)(x - r^2) = (x - i)(x + 1) = x^2 + (1 - i)x - i
    complex zero_poly[] = {{0, -1}, {1, -1}, {1, 0}, {0, 0}};
    complex zero_poly_eval[4];
    eval_from_poly(zero_poly_eval, zero_poly);
    print_poly("ZeroPoly eval", zero_poly_eval);

    // Create the evaluation (E * Z)(r^j) for each j. It has zeros at indices where we lost data.
    complex ez_eval[4];
    for (int j = 0; j < 4; j++) {
        ez_eval[j] = mul(data_eval[j], zero_poly_eval[j]);
    }
    print_poly("EZ eval", ez_eval);

    // Interpolate to get (E * Z)(x) = (D * Z)(x). The equality holds since, by construction, E and D agree both where
    // non-zero, and where zero due to Z.
    complex dz_poly[4];
    poly_from_eval(dz_poly, ez_eval);
    print_poly("EZ = DZ poly", dz_poly);

    // Scale ("shift") the polynomials so that we can divide them without hitting a zero in the zero poly.
    // Any scale factor is ok, as long as it is not one of the roots of unity (or zero).
    scale(dz_poly, two);
    scale(zero_poly, two);
    print_poly("DZ poly scaled", dz_poly);
    print_poly("ZeroPoly scaled", zero_poly);

    // Now we will divide the scaled polynomial (D * Z)(x / 2) by Z(x / 2), which will result in D(x / 2) - our
    // (scaled) original data. We do this via convolution: convert to evaluation form, divide pointwise, and convert
    // back to polynomial form.

    // Convert (D * Z)(x / 2) and Z(x / 2) to evaluation form,
    complex dz_scaled_eval[4];
    complex zero_scaled_eval[4];
    eval_from_poly(dz_scaled_eval, dz_poly);
    eval_from_poly(zero_scaled_eval, zero_poly);
    print_poly("DZ eval scaled", dz_scaled_eval);
    print_poly("Zero eval scaled", zero_scaled_eval);

    // Divide pointwise
    complex quotient_eval[4];
    for (int j = 0; j < 4; j++) {
        quotient_eval[j] = div(dz_scaled_eval[j], zero_scaled_eval[j]);
    }
    print_poly("Quotient eval", quotient_eval);

    // Convert back to polynomial form to get D(x / 2).
    complex recovered_poly[4];
    poly_from_eval(recovered_poly, quotient_eval);
    print_poly("Scaled recovered", recovered_poly);

    // Reverse our earlier scaling to recover D(x).
    unscale(recovered_poly, two);
    print_poly("Recovered values", recovered_poly);

    return 0;
}
