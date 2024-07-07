#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen
# add #include <algorithm> in vertex.cpp and CNTCell.cpp

cd dts_src
cl /O2 /EHsc /std:c++17 *.cpp
cl /FeDTS *.obj
move DTS.exe ..\.

cd ..
cd dts_convert
cl /O2 /EHsc /std:c++17 *.cpp
cl /FeCNV *.obj
move CNV.exe ..\.

cd ..
cd dts_generate
cl /O2 /EHsc /std:c++17 *.cpp
cl /FeGEN *.obj
move GEN.exe ..\.
cd ..