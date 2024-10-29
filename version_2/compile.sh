#!/bin/bash

# Script to compile DTS, CNV, and GEN with OpenMP if available
# Author: Weria Pezeshkian
# Niels Bohr International Academy, Niels Bohr Institute, University of Copenhagen

set -e  # Exit immediately on error

# Variables
SOURCE_DIR="dts_src"
CONVERT_DIR="dts_convert"
GENERATE_DIR="dts_generate"
CXX="g++-14"  # Set your specific compiler here
CXX_FLAGS="-O3 -std=c++11"
OPENMP_FLAG="-fopenmp"

# Print compiler version for confirmation
echo "Using compiler: $CXX"
$CXX --version

# Function to check for OpenMP support
check_openmp_support() {
    echo "Checking for OpenMP support with $CXX..."

    # Create a temporary file to test OpenMP support
    temp_file=$(mktemp /tmp/test_openmp.XXXXXX.cpp)
    echo '#include <omp.h>' > "$temp_file"
    echo "int main() { return 0; }" >> "$temp_file"
    
    # Try compiling the test file with OpenMP
    if $CXX $OPENMP_FLAG "$temp_file" -o /dev/null &>/dev/null; then
        echo "OpenMP support detected."
        rm "$temp_file"
        return 0
    else
        echo "OpenMP support not detected."
        rm "$temp_file"
        return 1
    fi
}

# Function to compile a module
compile_module() {
    local dir=$1
    local output_file=$2

    echo "Compiling in directory: $dir"
    cd "$dir" || { echo "Directory '$dir' not found"; exit 1; }

    # Compile based on OpenMP availability
    if check_openmp_support; then
        echo "Compiling $output_file with OpenMP support..."
        $CXX -c $CXX_FLAGS $OPENMP_FLAG *.cpp
        $CXX $OPENMP_FLAG -o "$output_file" *.o
    else
        echo "Compiling $output_file without OpenMP support..."
        $CXX -c $CXX_FLAGS *.cpp
        $CXX -o "$output_file" *.o
    fi

    # Move compiled output to parent directory
    mv "$output_file" ../
    echo "$output_file compiled and moved to parent directory."

    # Return to original directory
    cd ..
}

# Compile DTS, CNV, and GEN
compile_module "$SOURCE_DIR" "DTS"
compile_module "$CONVERT_DIR" "CNV"
compile_module "$GENERATE_DIR" "GEN"

# Completion message
echo "All modules compiled successfully. Executables DTS, CNV, and GEN are in the parent directory."
