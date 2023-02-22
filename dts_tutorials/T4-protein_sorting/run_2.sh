#!/bin/bash

rm top.top
echo "topol.q  10" > top.top
../../DTS -in input_2.dts -top top.top 
