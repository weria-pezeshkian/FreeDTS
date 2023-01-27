

#include <time.h>
#include "src/Analysis.h"



Analysis::Analysis(State *pState)
{
    m_TRJ = pState->GetTrajectory();
    
    
    
    std::cout<<m_TRJ.size()<<"size of trj \n";

    
    
    
    
}
Analysis::~Analysis()
{
}
