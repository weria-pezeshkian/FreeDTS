#!/bin/bash
../../GEN -box 50 50 100 -type flat -o topol.q
echo "topol.q  10" > top.top

../../DTS -in input.dts -top top.top 
