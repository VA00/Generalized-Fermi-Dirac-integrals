Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 9 Feb 2022

Reference arbitrary precision (WorkingPrecision->64) values for generalized relativistic Fermi-Dirac integrals (GFDI) computed by Mathematica 12.3
on 56-core  dual Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz, 376 GB memory. 

k=1/2

eta,theta = Full IEEE754 (double) test range for eta and theta in powers of 2, step 11.

-2 <= eta <= 8388608 = 2^23
0 <= theta <= 2^364 = 3.757668132438133e109

Any doubtful results were discarded. Timeout of 10 minutes, no memory constraints.
Mma code attached, refVALS.m.

To re-run use:

math <refVALS.m

Started: 09.02.2022 11:19
Running time: 11 hours  15 minutes  33 seconds
Finished: 09.02.2022 19:36

Input values: ?
After overflow check: ?
Actually computed: 24258

Files:

Mma_good_REF_finally_2022-02-09.m - Mma file with Rational numbedr inputs and 64-digit arbitrary precision output
refTBL_double_2022-02-09.bin - binary double precision dataset
refTBL_quad_2022-02-09.bin   - binary quadruple precision dataset


