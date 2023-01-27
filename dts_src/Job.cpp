

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
 Copyright (c) Weria Pezeshkian
 This class finds checks the name of the executable, Not a very important task for the current version
 */
#ifdef _OPENMP
# include <omp.h>
#endif
#include "Job.h"
#include "State.h"
#include "Nfunction.h"
#include "MC_Simulation.h"
#include "RNG.h"


Job::Job(std::vector <std::string> argument)
{
std::string Exe=argument.at(0);
    std::string ExeName = "DTS";
    Nfunction f;

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
#ifdef _OPENMP
State tem_S(argument);
Parallel_Tempering IPTdata = tem_S.m_Parallel_Tempering;
int Parallel_Tempering = IPTdata.State;
// if the state of Parallel_Tempering is off we move on to sim class, if not we do other things
if (Parallel_Tempering==0)
{
    if(tem_S.m_Integrator == "MC")
    {
        MC_Simulation SIM(&tem_S);
    }
}
else
{
            //  Parallel_Tempering  = on PT_steps  PT_minbeta    PT_maxbeta
            double PT_minbeta = IPTdata.PT_minbeta;          // from inputfile
            double PT_maxbeta = IPTdata.PT_maxbeta;            // from inputfile
            int PT_steps = IPTdata.PT_steps;             // from inputfile
            std::string gfile = tem_S.m_GeneralOutputFilename;
            std::string tsifile = (tem_S.m_TRJTSI).tsiFolder_name;
            RNG Random1(tem_S.m_Seed);
            int itime = tem_S.m_Initial_Step;
            int etime = tem_S.m_Final_Step;
            int num_threads=tem_S.m_Total_no_Threads;
            int exchange_step_length = (etime-itime)/PT_steps;
#if TEST_MODE == Enabled
std::cout<<"Parallel_Tempering "<<PT_steps<<" min beta "<<PT_minbeta<<" max beta "<<PT_maxbeta<<"\n";
#endif
            
           // int no_proc = omp_get_num_procs(); not sure about this
            //if(no_proc<num_threads)
            //num_threads = no_proc;
            omp_set_num_threads(num_threads);         // setting the number of thread
            std::vector<double> threads_energy;      // a container to store energy of each thread just (i) is the thread id
            std::vector<int>  tempid_thread_id;   // a map for temp_id to thread_id, temp = temperatures.at(temp_id)
            std::vector<double> betas;        // the 1/temprature of each temp_id
            for (int i=0;i<num_threads;i++)
            {
                threads_energy.push_back(0);
                betas.push_back(PT_minbeta + double(i)*(PT_maxbeta - PT_minbeta)/double(num_threads-1));
                tempid_thread_id.push_back(i);

            }
#pragma omp parallel if(Parallel_Tempering)
{             //firstprivate()
                State S(argument);
                int Thread_ID = omp_get_thread_num();
                int Thread_num = omp_get_num_threads();
    if(Thread_num!=num_threads)
    {
        std::cout<<" error ---> requested thread number is not provided \n";
        exit(0);
    }
                if(num_threads>1)
                S.m_Beta = PT_minbeta + double(Thread_ID)*(PT_maxbeta - PT_minbeta)/double(Thread_num-1);
                else
                {
                    std::cout<<" error--> Parallel_Tempering is ON while number of thread is 1 \n";
                    exit(0);
                }
    
                if(Thread_ID == Thread_num-1)
                    S.m_Targeted_State = true;
                else
                    S.m_Targeted_State = false;


#if TEST_MODE == Enabled
#pragma omp critical
std::cout<<"thread id "<<Thread_ID<<" total no threead "<<omp_get_num_threads()<<" no of assigned thread "<<m_Total_no_Threads<<"\n";
#endif
            //============================
            //========: Runing simulation
            for (int pt_step = 0; pt_step<exchange_step_length;pt_step++)
            {
                    S.m_Initial_Step = itime + pt_step*PT_steps;
                    S.m_Final_Step = itime + (pt_step+1)*PT_steps;
                    
                    // set the tempratur eof each state
                if(S.m_Integrator == "MC")
                {
                    MC_Simulation SIM(&S);
                    #pragma omp critical //(filling) not sure if it is needed
                    {
                        threads_energy[Thread_ID] = S.m_TotEnergy;   // get the latest energy of the systems
                    }
                    #pragma omp barrier
                    //change temprature
                    #pragma omp single
                    {
                        for (int c=0;c<betas.size()-1;c++)
                        {
                            double b1 = betas[c];
                            double b2 = betas[c+1];
                            int t1 = tempid_thread_id[c];
                            int t2 = tempid_thread_id[c+1];
                            double e1 = threads_energy[t1];
                            double e2 = threads_energy[t2];
                            double cra = (e2-e1)*(b2-b1);  // should be checked
                            double ran = Random1.UniformRNG(1.0); //should be checked
                            if(exp(cra)>ran)
                            {
                                // P = min(1,exp([E_i-Ej]*(1/Ti-1/Tj))) = min(1,exp([E_i-E_j]*(1/Ti-1/Tj)))
                                tempid_thread_id[c] = t2;
                                tempid_thread_id[c+1] = t1;
                            }
                        }
                    }
                    #pragma omp critical //(filling) not sure if it is needed
                    {
                        std::vector<int>::iterator it= find(tempid_thread_id.begin(), tempid_thread_id.end(), Thread_ID);
                        int temid_of_thread = std::distance(tempid_thread_id.begin(), it);
                         if(temid_of_thread<betas.size())
                         {
                             S.m_Beta = betas[temid_of_thread];
                             if(temid_of_thread == num_threads-1)
                                 S.m_Targeted_State = true;
                             else
                                 S.m_Targeted_State = false;
                         }
                        else
                        {
                            std::cout<<"error-->3627 should not happen\n";
                            exit(0);
                        }
                        
                    }
                    #pragma omp barrier
                    S.m_GeneralOutputFilename = gfile+ +"_temp_" + f.Int_to_String(S.m_Beta);
                    (S.m_TRJTSI).tsiFolder_name =  tsifile+"_temp_" +f.Int_to_String(S.m_Beta);
                }
                
            }
} // end pragma omp parallel
}  // end of Parallel_Tempering
#else   // When we do not use OPENMP
            {
                State S(argument);
                if(S.m_Integrator == "MC")
                {
                    MC_Simulation SIM(&S);
                }
                
            }
#endif
        }
        else
        std::cout<<" Error:  (c) unrecognized exacutable name --->"<<Exe<<" :( "<<std::endl;
    }


}
Job::~Job()
{
    
}

