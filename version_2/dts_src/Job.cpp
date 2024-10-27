/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#ifdef MPI_DETECTED
    #include <mpi.h> 
#endif
#include <vector>
#include <string>
#include "SimDef.h"
#include "Job.h"
#include "State.h"
#include "RNG.h"
#include "AbstractParallelReplicaRun.h"
#include "ParallelTempering.h"

/*
Description:
    This class handles task distribution and allows for execution.

Parameters:
    argument (std::vector<std::string>&): Vector containing input arguments.

*/
Job::Job(const std::vector<std::string> &argument) {
    // Extract executable name from the argument list
    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    
    // Check if the executable name matches the expected name
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }
#ifndef _OPENMP
    // Perform a normal simulation on a single CPU if OpenMP is not enabled
    State T_state(argument);
    T_state.Initialize();
    T_state.GetSimulation()->do_Simulation();
    
#else
//---> constract an State object
    State T_state(argument);
//---> get parallel tempering data in the input.dts file
    ParallelReplicaData    PRD = T_state.GetParallelReplicaData();
    
//---> here is one openmp is on but still want to perform one single simulation
    if (!(PRD.State)) {
        T_state.Initialize();
        T_state.GetSimulation()->do_Simulation();
    }
else { // run parallel tempering simulations
    AbstractParallelReplicaRun *pParallelReplicaRun;
        
    if (!(PRD.Type == ParallelTempering::GetDefaultReadName())){
            
        pParallelReplicaRun = new ParallelTempering(argument);
        if(pParallelReplicaRun->Initialize(PRD)){
            pParallelReplicaRun->Run();
        }
        else{
            std::cout<<"---> error: faild .... "<<"\n";
        }
    }
    else{
        std::cout<<"---> error: unknow type for "<<AbstractParallelReplicaRun::GetBaseDefaultReadName()<<"\n";
        exit(0);
    }
    }
#endif


}
Job::~Job() {
    
}




