#!/usr/bin/env bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_ana
g++ -c -O3 *.cpp
g++ -o ANA *.o
mv ANA ../
cd ..
