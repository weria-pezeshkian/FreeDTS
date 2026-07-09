/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#include <vector>
#include <string>
#include "SimDef.h"
#include "Job.h"
#include "State.h"
#include "RNG.h"


/*
Description:
    This class handles task distribution and allows for execution.

Parameters:
    argument (std::vector<std::string>&): Vector containing input arguments.

*/
//Just a little test
//Second test
Job::Job(const std::vector<std::string> &argument) {
    // Extract executable name from the argument list
    
    
    // Check if the executable name matches the expected name

    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }

    std::cout << "Running simulation on single CPU" << std::endl;
    State T_state(argument);
    T_state.Initialize();
    T_state.GetSimulation()->do_Simulation();
}
Job::~Job() {
    
}




