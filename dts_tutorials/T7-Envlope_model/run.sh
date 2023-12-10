#!/bin/csh


set ts2cg = '/Users/weria/Documents/TS2CG1.2/TS2CG1.2'
$ts2cg/PLM -TSfile  dts/TrjTSI/output10.tsi  -bilayerThickness 4.05 -rescalefactor 4.8 4.8  4.8   -Mashno 4 
$ts2cg/PCG  -dts point -str input.str -seed 512498  -Bondlength 0.2 -LLIB Martini3_all_lipids.LIB -Rcutoff 0.3 





  






