#!/bin/bash

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

cd src
g++ -c -O3 *.cpp
cd ..
g++ -c -O3 Analysis.cpp
g++ -o ANA Analysis.o src/*.o
rm src/*.o
rm *.o


