Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 9 Feb 2022

Reference values for generalized relativistic Fermi-Dirac integrals (GFDI) computed by Mathematica 12.3
on 56-core  dual Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz, 376 GB memory. 

k=0

eta,theta = Full IEEE754 (double) test range for eta and theta in powers of 2, step 11.

-256 <= eta <= 8388608
0 <= theta 2^353 = 1.8347988927920572e106

Any doubtful results were discarded. Timeout of 10 minutes, no memory constraints.
Mma code attached, refVALS.m.

To re-run use:

math <refVALS.m

Started: 08.02.2022 11:19
Running time: 18 hours  43 minutes  20 seconds
Finished: 09.02.2022 06:03

Input values: 35908
After overflow check: 35908
Actually computed: 7683	

Files:

Mma_good_REF_finally_2022-02-09.m - Mma file with Rational numbedr inputs and 64-digit arbitrary precision output
refTBL_double_2022-02-09.bin - binary double precision dataset
refTBL_quad_2022-02-09.bin   - binary quadruple precision dataset


