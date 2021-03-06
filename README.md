# Toy example of data recovery

The code implements, for didactic purposes, a very simple example of data recovery. It encodes a tiny case of Vitalik Buterin's method described in [this article](https://ethresear.ch/t/reed-solomon-erasure-code-recovery-in-n-log-2-n-time-with-ffts/3039?u=benjaminion).

See the accompanying write-up [here](https://hackmd.io/@benjaminion/data_recovery).

For something more capable, see my [`c-kzg` repo](https://github.com/benjaminion/c-kzg).

## Compile

Either

```
gcc main.c
```

or

```
clang main.c
```

Or choose your favourite C compiler.

## Run

```
./a.out
```

## Output

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
