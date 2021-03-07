# Toy example of data recovery

The code implements, for didactic purposes, a very simple example of data recovery. It exemplifies a tiny case of Vitalik Buterin's method described in [this article](https://ethresear.ch/t/reed-solomon-erasure-code-recovery-in-n-log-2-n-time-with-ffts/3039?u=benjaminion).

See the accompanying write-up [here](https://hackmd.io/@benjaminion/data_recovery). For something more capable, see my [`c-kzg` repo](https://github.com/benjaminion/c-kzg).

There are two versions:

1. A version that encodes the data and does the calculations over the complex integers. This uses _i_ as the primitive fourth root of unity, giving the roots for the Fourier transforms as [1, i, -1, -i]. Find it in _complex.c_.
2. A version that encodes the data and does the calculations over the integers modulo 17, a finite field. This uses 4 as the primitive fourth root of unity, giving the roots for the Fourier transforms as [1, 4, 16, 13]. Find it in _finite.c_.

## Compile

Either

```
cc complex.c
```

or

```
cc finite.c
```

Or choose your favourite C compiler.

## Run

```
./a.out
```

## Output

Note that, if you substitute `i = 4` in the output from _complex.c_, and make everything positive modulo 17, you end up with exactly the same output as from _finite.c_.

From _complex.c_

```
    Initial values: [5, 7, 0, 0]
      Data encoded: [12, 5 + 7i, -2, 5 - 7i]
 Data with missing: [12, 0, 0, 5 - 7i]
     ZeroPoly eval: [2 - 2i, 0, 0, -2 - 2i]
           EZ eval: [24 - 24i, 0, 0, -24 + 4i]
      EZ = DZ poly: [0 - 5i, 5 - 12i, 12 - 7i, 7]
    DZ poly scaled: [0 - 5i, 10 - 24i, 48 - 28i, 56]
   ZeroPoly scaled: [0 - i, 2 - 2i, 4, 0]
    DZ eval scaled: [114 - 57i, -24 - 23i, -18 - 9i, -72 + 69i]
  Zero eval scaled: [6 - 3i, -2 + i, 2 + i, -6 - 3i]
     Quotient eval: [19, 5 + 14i, -9, 5 - 14i]
  Scaled recovered: [5, 14, 0, 0]
  Recovered values: [5, 7, 0, 0]
```

From _finite.c_

```
    Initial values: [5, 7, 0, 0]
      Data encoded: [12, 16, 15, 11]
 Data with missing: [12, 0, 0, 11]
     ZeroPoly eval: [11, 0, 0, 7]
           EZ eval: [13, 0, 0, 9]
      EZ = DZ poly: [14, 8, 1, 7]
    DZ poly scaled: [14, 16, 4, 5]
   ZeroPoly scaled: [13, 11, 4, 0]
    DZ eval scaled: [5, 3, 14, 0]
  Zero eval scaled: [11, 2, 6, 16]
     Quotient eval: [2, 10, 8, 0]
  Scaled recovered: [5, 14, 0, 0]
  Recovered values: [5, 7, 0, 0]
  ```
