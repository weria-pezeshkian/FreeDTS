#include <iostream>
#include "SimDef.h"
#include "Job.h"
//To avoid having to include mpi.h, one can include a MACRO in a makefile that detects
//MPI...
#ifdef MPI_DETECTED
    #include <mpi.h>
#endif
/*
==============================================================================
 [FreeDTS] Analysis of Dynamically Triangulated Surfaces 
==============================================================================

Author: Weria Pezeshkian 2016-Now (Last Update 2024)
Contact: weria.pezeshkian@gmail.com

Description:
------------
FreeDTS is a computational tool designed for simulating biomembranes using dynamically triangulated surfaces. This code implements various algorithms and techniques to model membranes and their shapes.

Main Function:
--------------
This file contains the main function of the FreeDTS code. It initializes the simulation environment, reads input arguments from the command line, and passes them to the Job class responsible for executing the simulation tasks.

Usage:
------
To run FreeDTS, execute the compiled executable file with appropriate command-line arguments.
The executable should have the name of DTS, otherwise it does not run.
The name can be changed in the SimDef.h file
You can run "DTS -h" to get help about what kind of command-line arguments should be used.
Altenratively check the manual
 

Disclaimer:
-----------
This software is provided "as is," without warranty of any kind.
For more information and updates, please visit the project repository:
 https://github.com/weria-pezeshkian/FreeDTS

---------------------------------------------------------------------------
*/
int main(int argc, char* argv[]) {
    // Print program header
    
    // First, check if both _OPENMP and MPI are defined
    #if defined(_OPENMP) && defined(MPI_VERSION)
    #error "Both _OPENMP and MPI are defined and which have not been eimplemented. Exiting."
    // In the future, this might be the way to implement a hybrid OpenMP & MPI option.
    // Next, check if only _OPENMP is defined
    #elif defined(_OPENMP)
    // _OPENMP detected
    std::cout << "████████████████████████████████    FreeDTS    ██████████████████████████████" << std::endl;
    std::cout << "------ Analysis of Dynamically Triangulated Surfaces ------------" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------" << std::endl;
    std::cout << "------ OpenMP is supported and detected ----------------------------------------------" << std::endl;
    // Store command-line arguments in a vector
    std::vector<std::string> argument(argv, argv + argc);
    // Create a Job object and pass the arguments to it
    Job job(argument);

    #elif defined(MPI_VERSION)
    //MPI detected
    //Macro to check if MPI is detected 
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0){
    std::cout << "████████████████████████████████    FreeDTS    ██████████████████████████████" << std::endl;
    std::cout << "------ Analysis of Dynamically Triangulated Surfaces ------------" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------" << std::endl;
    std::cout << "------ MPI is supported and detected ----------------------------------------------" << std::endl;}
    std::vector<std::string> argument(argv, argv + argc);
    // Create a Job object and pass the arguments to it
    Job job(argument);

    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    #else
    std::cout << "████████████████████████████████    FreeDTS    ██████████████████████████████" << std::endl;
    std::cout << "------ Analysis of Dynamically Triangulated Surfaces ------------" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------" << std::endl;
    // Store command-line arguments in a vector
    std::vector<std::string> argument(argv, argv + argc);
    // Create a Job object and pass the arguments to it
    Job job(argument);
    #endif
    
    

    

    return 0;
}
