/*
 Simulation of biomembranes using dynamically triangulated surfaces:
 by Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 
 
 This is the main class of the code: This class receives the input arguments and pass it on to the job class.
 */
#include "SimDef.h"
#include "Job.h"

int main(int argc, char* argv[])
{

    std::cout << "----------------- "<<S_Name<<" ---------------- "<<std::endl;
    std::cout << "------ Simulation of Biomembranes Using Dynamically Triangulated Surfaces ------------ "<<std::endl;
    std::cout << "---------------------------- Weria Pezeshkian -------------------------------------- "<<std::endl;

    std::vector <std::string> argument;
    std::string Temprory;
           for (long i=0;i<argc;i++)
           {
               Temprory.assign(argv[i]);
               argument.push_back(Temprory);
           }
           Job job(argument);


    return 0;
}
