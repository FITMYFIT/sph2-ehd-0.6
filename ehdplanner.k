*KEYWORD
$
$ NODES
$
*RUN_MODE 0-Implicit 1-Explicit 
1
$
*START_STEP 0-from beginning; else-from specified step 
0
$
*PART 1PID 2PTTYPEID 3MID 4EOSID 5C0 6HdivDp 7VISK 8VISN 9EKAPPA 10 EEPSILON
1,2,1,2,1,1.5,1.0e-3,0.27,0,3
2,1,1,2,1,1.5,1.0e-3,0.27,0,3
3,1,1,2,1,1.5,1.0e-3,0.27,0,1
4,2,1,2,1,1.5,1.0e-3,0.27,0,1
$
*CONTROL_SPH1 1SPHST 2SPHPV 3SPHAV 4SPHAS 5SPHAD 6XSPH-Epsilon 7Period-Boundary 8Damp-T 9Particle-Shift-Coef
0,0,2,0,0,0.0,0,0,0
$
*CONTROL_SPH2 1Sound-Speed 2DeltaT 3FinalTime 4OutputSteps 
1.0,1.0E-6,1,1
$
*CELL_NUM 1CELLNUMX,2CELLNUMY
50,50
$
*CALCULADimGreyTION_RANGE 1Xmin 2Xmax,3Ymin 4Ymax
0.0 0.0,0.0 0.0
$
*SPH_AV 1AVAlpha 2AVBeta 3AVEta
0.02,0.01,0.1
$
*SPH_AS 1ASEpsilon1 2ASEpsilon2 3ASDeltaD
0.2,0.3,2.0e-5
$
*DENSITY_RENORMSTEPS 0-No Renormlize else-Renormlized Steps
0
$
*PERIOD_BOUNDARY 1xmin 2xmax
-2.0e-4 2.0e-4
$
*EOS_IDEAL_GAS 1EOSID 2Cp 3Mw 4R
1,14.545,2.0,8.31
$
*EOS_WEAKLY_COMPRESS 1EOSID 2P0 3GAMMA 4EXTP
2,2000,7.0,0.0
$
*SURFACE_TENSION 1PID 2STSigma 3STHvsh(Must be Same as Part Num)
1,0.0727,1.5
2,0.0727,1.5
3,0.0727,1.5
$
*LOAD_BODYFORCE 1PID 2fx 3fy
1,0.0,0.0
2,0.0,0.0
3,0.0,0.0
$
*INITIAL_SPH 1PID 2u 3v 4rho 5T
1,0.0,0.0,1000.0,0.0
2,0.0,0.0,1000.0,0.0
3,0.0,0.0,1000.0,0.0
$
*END