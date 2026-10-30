#ifndef JOB_FREEDTS_H_INCLUDED
#define JOB_FREEDTS_H_INCLUDED

#include "SimDef.h"
#include <vector>
#include <string>

class Job {
public:
    // Constructor
    Job(const std::vector<std::string>& argument);
    
    // Destructor
    ~Job();
    
private:
    // Private member variables and functions can be added here
};

#endif // JOB_H_INCLUDED
