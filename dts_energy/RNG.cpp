


#include "RNG.h"
#include "Nfunction.h"
RNG::RNG(int seed)
{
m_Seed=seed;
    Nfunction f;
    
#if RNGTYPE == UNIFROMTYPE1
    // this part have not been compeleted. not recommended to be used
    std::default_random_engine generator (m_Seed);
    m_generator = generator;
    std::uniform_real_distribution<double> Rdistribution (0.0,1.0);
    m_Rdistribution = Rdistribution;
   // std::uniform_int_distribution<int> distribution(1,6);
    std::string A=" RNG type 1 is used";
    f.Write_One_LogMessage(A);
#elif RNGTYPE == UNIFROMTYPE0
    srand (m_Seed);
    std::string A=" RNG type 0 is used";
    f.Write_One_LogMessage(A);
#endif

}


RNG::~RNG()
{
    
}

 
double RNG::UniformRNG(double A)
{
double rng=0.0;

#if RNGTYPE == UNIFROMTYPE1
   rng = A*(m_distribution(m_generator));
#elif RNGTYPE == UNIFROMTYPE0
    rng=A*(double(rand()%2000000)/2000000.0);
#endif



return rng;
}
int RNG::IntRNG(int a)
{
    int rng=0;
    int zz =0;
    
#if RNGTYPE == UNIFROMTYPE1
    zz = int(a*(m_distribution(m_generator)));
#elif RNGTYPE == UNIFROMTYPE0
    zz=rand()%a;
#endif
    
    
    if(zz==a)
    {
      rng = a-1; 
    Nfunction f;
    std::string B="Error 2. Random number should not happend";
    f.Write_One_LogMessage(B);   
    std::cout<<B<< "\n";
    }
    else
    {rng = zz;}
    
    return rng;
}


int RNG::BinRNG()
{
    int rng=0.0;
    
    double zz = UniformRNG(1)-0.5;
    if(zz>0)
    {rng = 1;}
    else
    {rng = -1;}
    
    return rng;
}





