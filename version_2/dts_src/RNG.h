#if !defined(AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234_INCLUDED_)
#define AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234_INCLUDED_

#include <random>
#include <ctime>
#include "SimDef.h"
/*
    RNG.h
    Random Number Generator (RNG) class for generating pseudo-random numbers.
    The RNG class provides methods for generating uniformly distributed random numbers,
 Note:
 - The behavior of the RNG class may vary depending on the preprocessor macros defined in the SimDef.h file.
 - This class uses either the C standard library function `rand()` for random number generation or
   the `<random>` library for more advanced random number generation, based on the defined macro RNGTYPE.
 - When using RNGTYPE UNIFROMTYPE1, the `<random>` library is utilized, providing better quality random numbers.
 
    Usage:
    - Construct an instance of RNG with a seed value to initialize the random number generator.
    - Call the UniformRNG method to generate a uniformly distributed random number between 0 and 1,
      or specify a maximum value to scale the result.
    - Call the IntRNG method to generate a random integer within a specified range.
    - Call the BinRNG method to generate a binary random value.

    Example:
    // Initialize RNG with a seed value
        RNG rng(123);

    // Generate a uniformly distributed random number between 0 and 1
    double randomValue = rng.UniformRNG();

    // Generate a random integer between 0 and 10
    int randomInt = rng.IntRNG(10);

    // Generate a binary random value
    bool binaryValue = rng.BinRNG();
*/
class RNG {
public:
    // Constructor: Initialize RNG with a seed
    explicit RNG(int seed);

    // Generates a uniform random number in the range [0, A)
    double UniformRNG(double A);
	inline int GetSeed()        {return m_seed;}
    // Generates a random integer in the range [0, a)
    int IntRNG(int a);

    // Generates a binary random outcome (true/false)
    bool BinRNG();

private:
    int m_seed;  // Seed for the random number generator

    // C++11 random number generator and distribution (used if RNGTYPE == UNIFORMTYPE1)
#if RNGTYPE == UNIFORMTYPE1
    std::default_random_engine m_generator;  // Random number generator (C++11)
    std::uniform_real_distribution<double> m_uniform_dist;  // Uniform real distribution [0, 1]
#endif

    // Note: No additional member variables are required for RNGTYPE == UNIFORMTYPE0 since it uses rand()/srand().
};

#endif
