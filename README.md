This is a PoC implementation of CTIDH-1024 and CTIDH-2048, with the original and our refined group action evaluation. 

It is not meant for direct comparison with the original CTIDH implementation. When integrated into the CTIDH, these minor refinements may only improve the performance by **2~4%**. 

But **when allowing parallel, the refined version can be significantly faster than the original**. See remarks below. _Not interesting, but it can be useful._

## Dependencies

- SageMath (v10.1+ for full tests; v9.5+ should suffice for benchmarking)
- Additional dependencies: `pip install -r requirements.txt`

*Key dependencies include:*
- Sage's finite field arithmetic (the `Primefield` wrapper only for counting operations)
- Elliptic curve and isogeny computation for testing / verification of correctness.


## Benchmarking
1. git checkout first:
   - For "original CTIDH", `git checkout source-code-CTIDH`.
   - For our refined CTIDH, `git checkout improved-CTIDH`
2. Then run `python3 -m unittest tests.test_csidh.TestCSIDH.test_group_action` with sagemath's python.

## Comparison with Original CTIDH Implementation
### Algorithmic Differences
This implementation differs from original CTIDH in:
1. Uses only traditional Vélu for isogeny computations (results in higher computational cost)
2. (Virtually) use Montgomery ladder for modular exponentiations (without window/NAF techniques)


### Group Action Evaluation Optimizations (branch `improved-CTIDH`)
1. **Single random point sampling:** Samples Elligator points just once per inner loop to reduce sampling cost
2. **Optimized batch selection with adjusted key spaces to minimize xEVAL costs:**
   - Constrains batch bounds for last few batches (sum of absolute values ≤ max batch bound)
   - Similarly constrains preceding batch bounds
   - Processes only one final batch per inner loop (saves xEVALs for mentioned batches)
3. **Redundant scalar multiplication eliminated:** Removes scalar multiplications of the unused point when inner loop has single batch

**NOTE:**  Key spaces were manually adjusted from original CTIDH, **with some decrease of the sizes**. Key space sizes (in bits):

|      | `source-code-CTIDH` | `improved CTIDH` |
| ---- | -------------- | ---------------- |
| 1024 | 256.1          | **244.9**        |
| 2048 | 220.0          | **218.5**        |


<!-- ## Some Further remarks
If we allow parallelization, our refinements would be more valuable:
1. Samping point only once makes it possible to **compute the scalar multiplications of two random points in parallel**. These scalar mults are expensive.
2. Although the cost of random points sampling is only about 5% of the original CTIDH group action (without parallel), and original CTIDH can save a xEVAL for an extra isogeny. But it could be useful, because:
   1. The two xEVAL of an isogeny can be done in parallel, therefore saving a single xEVAL is not attractive
   2. While all other algorithms (i.e. scalar mult, isogeny) are highly parallelizable (at least, in theory), Elligator, whose cost mainly originated from Legendre symbol computation is not. -->

## Some Further remarks
1. Parallel Scalar Multiplication: **Sampling points just once enables parallel computation of two random points' scalar multiplications** - these are expensive.
2. Optimization Trade-offs Analysis:
While random point sampling constitutes only ~5% of original CTIDH's group action cost (serial execution),
and original CTIDH saves one xEVAL per extra isogeny. Our refinement remains valuable because: 
   - Parallel processing handles both xEVALs simultaneously, making single xEVAL savings less significant
   - **Other operations (such as scalar multiplication and isogeny computation) are highly parallelizable (in theory)**. Therefore Elligator's performance (dominated by Legendre symbol computation) remains serial-bound

## License
This repo is licensed under the GPL v3 (see LICENSE file).

*Note: Contains code derived from the [sibc](https://github.com/JJChiDguez/sibc) library.*
