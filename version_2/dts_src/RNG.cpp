
#include "RNG.h"
// Constructor
RNG::RNG(int seed) : m_seed(seed) {
    // Initialize the generator based on the RNG type
    #if RNGTYPE == UNIFORMTYPE1
        m_generator.seed(m_seed);  // Use the seed to initialize the generator
        m_uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);  // Uniform distribution between 0 and 1
    #elif RNGTYPE == UNIFORMTYPE0
        srand(m_seed);  // Seed rand() with the given seed
    #endif
}

// Uniform RNG generator for doubles in range [0, A)
double RNG::UniformRNG(double A) {
    #if RNGTYPE == UNIFORMTYPE1
        return A * m_uniform_dist(m_generator);  // Use C++11 uniform distribution
    #elif RNGTYPE == UNIFORMTYPE0
        return A * (static_cast<double>(rand()) / RAND_MAX);  // Use rand() scaled to [0, A)
    #endif
}

// Uniform RNG generator for integers in range [0, a)
int RNG::IntRNG(int a) {
    #if RNGTYPE == UNIFORMTYPE1
        std::uniform_int_distribution<int> int_dist(0, a - 1);
        return int_dist(m_generator);  // Use C++11 uniform int distribution
    #elif RNGTYPE == UNIFORMTYPE0
        return rand() % a;  // Use rand() to generate integer in [0, a)
    #endif
}

// Bernoulli RNG for a binary outcome (true/false)
bool RNG::BinRNG() {
    // Generate a uniform random number in range [0, 1] and return true if greater than 0.5
    return UniformRNG(1.0) > 0.5;
}
