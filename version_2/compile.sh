#!/bin/bash

set -e

# Weria Pezeshkian
# Niels Bohr International Academy
# Niels Bohr Institute
# University of Copenhagen

(
  cd dts_src
  g++ -c -O3 *.cpp
  g++ -o DTS *.o
)

mv dts_src/DTS .

(
  cd ../version_1/dts_convert
  g++ -c -O3 *.cpp
  g++ -o CNV *.o
)

mv ../version_1/dts_convert/CNV .

(
  cd ../version_1/dts_generate
  g++ -c -O3 *.cpp
  g++ -o GEN *.o
)

mv ../version_1/dts_generate/GEN .
