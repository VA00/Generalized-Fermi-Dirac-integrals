Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 14 Feb 2022

Reference arbitrary precision (WorkingPrecision->64) values for generalized relativistic Fermi-Dirac integrals (GFDI) computed by Mathematica 12.3
on 56-core  dual Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz, 376 GB memory. 

k=1/2

eta,theta = Full IEEE754 (double) test range for eta and theta in powers of 2, step 3.

-2 <= eta <= 8388608 = 2^21
0 <= theta <= 2^351 = 4.586997231980143e105

Any doubtful results were discarded. Timeout of 10 minutes, no memory constraints.
Mma code attached, refVALS.m.

To re-run use:

math <refVALS.m

Started: 15.02.2022 06:31
Running time: 2 days 7 hours  23 minutes
Finished: 17.02.2022 13:55

Input values: 241452
After overflow check: 9614
Actually computed correctly: 6411

Files:

Mma_good_REF_finally_2022-02-17.m - Mma file with Rational numbedr inputs and 64-digit arbitrary precision output
refTBL_double_2022-02-17.bin - binary double precision dataset
refTBL_quad_2022-02-17.bin   - binary quadruple precision dataset


