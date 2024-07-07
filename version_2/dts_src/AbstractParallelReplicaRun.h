#ifndef AbstractParallelReplicaRun_FREEDTS_H_INCLUDED
#define AbstractParallelReplicaRun_FREEDTS_H_INCLUDED

#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class AbstractParallelReplicaRun {
public:
    // Constructor
    AbstractParallelReplicaRun(){
        
    }
    // Destructor
    virtual ~AbstractParallelReplicaRun(){
        
    }
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "ParallelReplicaRun";}
    
    virtual bool Initialize(ParallelReplicaData PRD) = 0;
    virtual bool Run() = 0;
    
private:

};


#endif // AbstractParallelReplicaRun_H_INCLUDED
