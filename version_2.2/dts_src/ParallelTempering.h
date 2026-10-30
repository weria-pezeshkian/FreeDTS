#ifndef ParallelTempering_FREEDTS_H_INCLUDED
#define ParallelTempering_FREEDTS_H_INCLUDED

#include "AbstractParallelReplicaRun.h"
#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class ParallelTempering : public AbstractParallelReplicaRun {
public:
    // Constructor
    ParallelTempering(std::vector<std::string> Argument);
    
    // Destructor
    ~ParallelTempering();
    inline std::string GetDerivedDefaultReadName()  {return "ParallelTempering";}
    inline static std::string GetDefaultReadName()  {return "ParallelTempering";}
    std::string CurrentState();
    
    bool Initialize(ParallelReplicaData PRD);
    bool Run();
    
private:
    // Private member variables and functions can be added here
    int m_Rate;
    int m_Bins;
    double m_minBeta;
    double m_maxBeta;
    std::vector<std::string> m_Argument;
    
};

#endif // ParallelTempering_H_INCLUDED
