#!/bin/bash
../../../GEN -box 100 100 100 -type flat -o topol.q
echo "topol.q  10" > top.top

../../../DTS -in input_2.dts -top top.top
