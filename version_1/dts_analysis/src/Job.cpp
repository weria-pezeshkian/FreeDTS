

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class finds checks the name of the executable, Not a very important task for the current version
 */
#include "Job.h"
#include "State.h"

Job::Job(std::vector <std::string> argument)
{
std::string Exe=argument.at(0);
    std::string ExeName = "ANA";
 

	if (Exe.size()<3)
	{
        std::cout<<" Error:  (a) unrecognized exacutable name --->"<<Exe<<" :( "<<std::endl;
    }
    else if (Exe.size()>3 && Exe.at(Exe.size()-4)!='/')
    {
        std::cout<<" Error:  (b) unrecognized exacutable name --->"<<Exe<<" :( "<<std::endl;
    }
    else
    {
        char L1 = Exe.at(Exe.size()-1);
        char L2 = Exe.at(Exe.size()-2);
        char L3 = Exe.at(Exe.size()-3);
        if(L1==ExeName.at(2) && L2==ExeName.at(1) && L3==ExeName.at(0))
        {
            State S(argument);
        }
        else
        std::cout<<" Error:  (c) unrecognized exacutable name --->"<<Exe<<" :( "<<std::endl;
    }


}
Job::~Job()
{
    
}

