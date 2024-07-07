#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd dts_src
g++ -c -O3 *.cpp
g++ -o DTS *.o
rm *.o
mv DTS ../
cd ..
cd dts_convert
g++ -c -O3 *.cpp
g++ -o CNV *.o
rm *.o
mv CNV ../
cd ..

cd dts_generate
g++ -c -O3 *.cpp
g++ -o GEN *.o
rm *.o
mv GEN ../
cd ..



