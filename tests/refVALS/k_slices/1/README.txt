Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 12 Feb 2022

Reference arbitrary precision (WorkingPrecision->64) values for generalized relativistic Fermi-Dirac integrals (GFDI) computed by Mathematica 12.3
on 56-core  dual Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz, 376 GB memory. 

k=1

eta,theta = Full IEEE754 (double) test range for eta and theta in powers of 2, step delta=6.

-16 <= eta <= 4194304=2^22
0 <= theta 2^370 = 2.4049076047604052e111

Any doubtful results were discarded. Timeout of 10 minutes, no memory constraints.
Mma code attached, refVALS.m.

To re-run use:

math <refVALS.m

Started: 10.02.2022 19:04
Running time: 1 day 11 hours 52 minutes 43.217627 seconds
Finished: 12.02.2022 06:58

Input values:  177674
After overflow check: 119383
Actually computed: 80394	

Files:

Mma_good_REF_finally_2022-02-12.m - Mma file with Rational numbedr inputs and 64-digit arbitrary precision output
refTBL_double_2022-02-12.bin - binary double precision dataset
refTBL_quad_2022-02-12.bin   - binary quadruple precision dataset


