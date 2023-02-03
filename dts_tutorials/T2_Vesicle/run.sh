#!/bin/bash

../../GEN -box 50 50 50 -type tetrahedron -N 20 -o topol.q

echo "topol.q 22" > top.top

../../DTS -in input.dts -top top.top
