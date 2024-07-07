#if !defined(AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234__INCLUDED_)
#define AFX_RNG_H_9G4B21B8_C13Q_5648_BF23_124095086234__INCLUDED_
#include <random>
#include "SimDef.h"

class RNG
{
public:
    
	RNG(int);
	 ~RNG();



    


public:
    
       // inline const std::string GetFileName()    const {return m_FileName;}

private:
    double m_Seed;




private:
    
#if RNGTYPE == UNIFROMTYPE1
    std::default_random_engine m_generator;
    std::uniform_real_distribution<double> m_Rdistribution;
#endif

    
public:
    double UniformRNG(double );
    int BinRNG();
    int IntRNG(int );
    



    





};


#endif
