#include <iostream>
#include "SimDef.h"
#include "Job.h"
/*
==============================================================================
 [FreeDTS] Simulation of Biomembranes Using Dynamically Triangulated Surfaces 
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
    std::cout << "████████████████████████████████    FreeDTS    ██████████████████████████████" << std::endl;
    std::cout << "------ Simulation of Biomembranes Using Dynamically Triangulated Surfaces ------------" << std::endl;
    std::cout << "--------------------------------------------------------------------------------------" << std::endl;

    // Store command-line arguments in a vector
    std::vector<std::string> argument(argv, argv + argc);

    // Create a Job object and pass the arguments to it
    Job job(argument);

    return 0;
}
