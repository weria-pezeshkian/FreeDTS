

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
 Copyright (c) Weria Pezeshkian
 This class finds checks the name of the executable, Not a very important task for the current version
 */

#include "Job.h"
#include "State.h"
#include "Nfunction.h"
#include "Analyze_Energy.h"
#include "RNG.h"


Job::Job(std::vector <std::string> argument)
{
std::string Exe=argument.at(0);
    Nfunction f;



                State S(argument);
                Analyze_Energy ANA(&S);

}
Job::~Job()
{
    
}

