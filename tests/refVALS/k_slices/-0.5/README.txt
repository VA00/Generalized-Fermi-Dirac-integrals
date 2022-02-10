Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 9 Feb 2022

Reference arbitrary precision (WorkingPrecision->64) values for generalized relativistic Fermi-Dirac integrals (GFDI) computed by Mathematica 12.3
on 56-core  dual Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz, 376 GB memory. 

k=-1/2

eta,theta = Full IEEE754 (double) test range for eta and theta in powers of 2, step 11.

-2 <= eta <= 8388608 = 2^23
0 <= theta <= 2^353 = 1.8347988927920572e106

Any doubtful results were discarded. Timeout of 10 minutes, no memory constraints.
Mma code attached, refVALS.m.

To re-run use:

math <refVALS.m

Started: 09.02.2022 20:52
Running time: 21 hours 15 minutes 50 seconds
Finished: 10.02.2022 18:09

Input values: ?
After overflow check: ?
Actually computed: 18946

Files:

Mma_good_REF_finally_2022-02-10.m - Mma file with Rational numbedr inputs and 64-digit arbitrary precision output
refTBL_double_2022-02-10.bin - binary double precision dataset
refTBL_quad_2022-02-10.bin   - binary quadruple precision dataset


