#include "RNG.h"

RNG::RNG(int seed) : m_Seed(seed) {
}
void RNG::Initialize(){
    // Initialize random number generator based on RNGTYPE
#if RNGTYPE == UNIFROMTYPE1
    std::default_random_engine generator(m_Seed);
    m_generator = generator;
#elif RNGTYPE == UNIFROMTYPE0
    srand(m_Seed);
#endif
}
double RNG::UniformRNG(double A) {
#if RNGTYPE == UNIFROMTYPE1
    return A * (m_distribution(m_generator));
#elif RNGTYPE == UNIFROMTYPE0
    return A * (double(rand() % 2000000) / 2000000.0);
#endif
}

int RNG::IntRNG(int a) {
#if RNGTYPE == UNIFROMTYPE1
    return int(a * (m_distribution(m_generator)));
#elif RNGTYPE == UNIFROMTYPE0
    return rand() % a;
#endif
}

bool RNG::BinRNG() {
    // Generate a random number between 0 and 1, and return true if it's greater than 0.5
    return UniformRNG(1) > 0.5;
}
