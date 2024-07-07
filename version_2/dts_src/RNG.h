#if !defined(AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234_INCLUDED_)
#define AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234_INCLUDED_

#include <random>
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
    RNG(int seed);
    ~RNG(){
        
    }
    void Initialize();

    double UniformRNG(double A);
    bool BinRNG();
    int IntRNG(int a);

private:
    double m_Seed;

#if RNGTYPE == UNIFROMTYPE1
    std::default_random_engine m_generator;
    std::uniform_real_distribution<double> m_Rdistribution;
#endif
};

#endif
