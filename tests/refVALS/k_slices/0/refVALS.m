(* Andrzej Odrzywołek, andrzej.odrzywolek@uj.edu.pl, 9 Feb 2022 *)
(* Reference value generator for generalized relativistic Fermi-Dirac integrals (GFDI) *)
$IgnoreEOF=1;
$MaxExtraPrecision = 1024;
(* Main function *)
F[k_, \[Eta]_, \[Theta]_, m_, n_] := Module[{integrand, f, x},
  f = Function[{z, kk, \[Eta]\[Eta], \[Theta]\[Theta]}, 
    z^kk Sqrt[
      1 + \[Theta]\[Theta] z/2] LogisticSigmoid[\[Eta]\[Eta] - z]];
  integrand = (Derivative[0, 0, m, n][f])[x, k, \[Eta], \[Theta]];
  TimeConstrained[If[\[Eta] > 4,
    Check[
      NIntegrate[Evaluate@integrand, {x, 0, \[Eta]}, 
       Method -> "GlobalAdaptive", MinRecursion -> 13, 
       MaxRecursion -> 64, WorkingPrecision -> 64], Infinity] +
     Check[
      NIntegrate[Evaluate@integrand, {x, \[Eta], Infinity}, 
       Method -> "GlobalAdaptive", MinRecursion -> 11, 
       MaxRecursion -> 29, WorkingPrecision -> 64], Infinity],
    Check[
     NIntegrate[Evaluate@integrand, {x, 0, Infinity}, 
      WorkingPrecision -> 64, Method -> "DoubleExponential", 
      MaxRecursion -> 17], Infinity]
    ], 600 (* 10 minutes and give UP! *), Infinity]
  ];
(* MAIN reference values generator. Beware of Mathematica BUGS !!!!! *)
(* kTBL = Join[{-1/2, 0}, Table[2^i, {i, -1, 6}]];*)
delta=11;
kTBL={0};
\[Theta]TBL = Join[{0}, Table[2^i, {i, -1022, 1024, delta}]];
(* \[Theta]TBL = {1}; *)
\[Eta]TBL = 
  Join[Table[-2^i, {i, 8, -1022, -delta}], {0}, 
   Table[2^i, {i, -1022, 26, delta}]];
(*\[Eta]TBL=SetPrecision[Import["H:\\Fermi-Dirac\\Generalized-Fermi-Dirac-integrals\\tests\\continuity.dat","TSV"][[1;;,-2]],64];*)
inputTBLall=Flatten[Table[{k,\[Eta],\[Theta]},{k,kTBL},{\[Eta],\[Eta]TBL},{\[Theta],\[Theta]TBL}],2];
inputTBLall//Length
(* Exclude overflof/underflow test inputs *)
overflowPositiveEta[k_,\[Eta]_,\[Theta]_]:=And[\[Eta]>143,\[Eta]^(1+k)/(3/2+k)Sqrt[1+\[Eta] \[Theta]/2]>2^1024];
overflowNegativeEta[k_,\[Eta]_,\[Theta]_]:=And[\[Eta]<0,Exp[\[Eta]]Max[Gamma[1+k],Gamma[3/2+k]*Sqrt[\[Theta]/2]]>2^1024];
underflowNegativeEta[k_,\[Eta]_,\[Theta]_]:=And[\[Eta]<0,\[Eta]+LogGamma[1+k]+Log[1+Gamma[3/2+k]/Gamma[1+k]Sqrt[\[Theta]/2]]<Log[2^-1022]];
inputTBL=Select[inputTBLall,!Or[underflowNegativeEta@@#,overflowNegativeEta@@#,overflowPositiveEta@@#]&];
Length[inputTBL]
timerStart=AbsoluteTime[];
mma3=ParallelTable[Table[If[m+n<=3,F[inputTBL[[ii,1]],inputTBL[[ii,2]],inputTBL[[ii,3]],m,n],Nothing],{m,0,3},{n,0,3}]//Flatten,{ii,1,Length[inputTBL]}];
time3=AbsoluteTime[]-timerStart
refTBLtmp=Join[Transpose[inputTBL],Transpose[mma3]]//Transpose;
refTBL=Select[refTBLtmp, NumericQ[Total[#]] &];
Export["refTBL_double_"<>DateString["ISODate"]<>".bin",refTBL,"Real64"]
Export["refTBL_quad_"<>DateString["ISODate"]<>".bin",refTBL,"Real128"]
Export["Mma_good_REF_finally_"<>DateString["ISODate"]<>".m",refTBL]
Print["Total reference values:\t",Length[refTBL]];
Print["Total time:\t",time3];
Print[DateString[]];


