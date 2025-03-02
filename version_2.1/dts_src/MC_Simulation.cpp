

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include <ctime>
#include <iostream>
#include <string.h>
#include "MC_Simulation.h"
#include "State.h"
#include "SimDef.h"
/*
 List of skipped function due to lack of clarity based on current state of the code
 They need to be finished before calling this a new version.
 
 1) void TimeSeriesLogInformation::WriteStateInfo(){
 
 */






/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 MC simulation class, runs mc simulation if it is defined in the input file.
 */
MC_Simulation::MC_Simulation(State *pState) : m_pState(pState) {

}
MC_Simulation::~MC_Simulation(){
    
}
void MC_Simulation::Initialize(){
    
    return;
}
bool MC_Simulation::do_Simulation(){
#if DEBUG_MODE == Enabled
    std::cout<<" do_Simulation function is starting  \n";
#endif
    
//---> Voxelize the mesh for the first time; this should be done before any calculation
    m_pState->GetVoxelization()->Voxelize(m_pState->GetMesh()->GetActiveV());
    
#if DEBUG_MODE == Enabled
    std::cout<<" system has been voxelaized  \n";
#endif

//----> checking if the mesh is good, within the bond of the simulation type. For here, it should be within
        //CheckMesh();
 
//--- before simualtion lets have a frame of the initial system
        m_pState->GetVisualization()->WriteAFrame(0);
   // time_t startTime;
   // time(&startTime);
#if DEBUG_MODE == Enabled
    std::cout<<" We have reached simulation run loop!  \n";
#endif
    std::clock_t start = std::clock();
    
#ifdef _OPENMP
    double startwall_time = omp_get_wtime();
#endif
// Your OpenMP parallel code here

    if(!m_pState->GetMesh()->CheckMesh(m_MinLength2, m_MaxLength2, m_MinAngle, m_pState->GetVoxelization())){
        
        std::cout << "---> error: The  mesh quality is insufficient for running a simulation.\n";
        exit(0);
    }

    
    

    std::cout<<"------>   Simulation will be performed from "<<m_Initial_Step<<" to "<<m_Final_Step<<" steps\n";
for (int step = m_Initial_Step; step <= m_Final_Step; step++){
        
//----> write files
        //--- write visulaization frame
        m_pState->GetVisualization()->WriteAFrame(step);
        //--- write non-binary trejectory e.g., tsi, tsg
        m_pState->GetNonbinaryTrajectory()->WriteAFrame(step);
        //--- write binary trejectory e.g., bts
        m_pState->GetBinaryTrajectory()->WriteAFrame(step);
        //--- write into time seri file, e.g., energy, volume ...
        m_pState->GetTimeSeriesDataOutput()->WriteTimeSeriesDataOutput(step);
        //--- write check point for the state
        m_pState->GetRestart()->UpdateRestartState(step, m_pState->GetVertexPositionUpdate()->GetDR(), m_pState->GetDynamicBox()->GetDR());
    
//---> centering the simulation box
    if(m_CenteringFrequently != 0 && step%m_CenteringFrequently == 0){
        m_pState->GetMesh()->CenterMesh();
        m_pState->GetVoxelization()->ReassignMembersToVoxels(m_pState->GetMesh()->GetActiveV());
    } // [ if(GetBoxCentering()!=0 && step%GetBoxCentering()==0)]

//---> Run standard Integrators
        //--- run the vertex position update
        m_pState->GetVertexPositionUpdate()->EvolveOneStep(step); // we may need the final step as well to check if the update of move size should be done
        //--- run the link flip update
        m_pState->GetAlexanderMove()->EvolveOneStep(step);
        //--- run the inclusion update
        m_pState->GetInclusionPoseUpdate()->EvolveOneStep(step);
        //--- run vector fields
        m_pState->GetVectorFieldsRotationUpdate()->EvolveOneStep(step);

//----> Run the supplementary integrators
       //--- update the box side
         m_pState->GetDynamicBox()->ChangeBoxSize(step); // we may need the final step as well to check if the update of move size should be done
        //--- update edge of mesh open edge
         m_pState->GetOpenEdgeEvolution()->Move(step);
        //--- update the mesh topology
        m_pState->GetDynamicTopology()->MCMove(step);
        //---- convert inclusions
        m_pState->GetInclusionConversion()->Exchange(step);
    
        //---- NonequilibriumCommands
        m_pState->GetNonequilibriumCommands()->Run(step);

//----> print info about the simulation, e.g., rate,
   // time_t currentTime;
   // time(&currentTime);
    if(!CheckMesh(step)){
        std::cout<<"---> error, the mesh does not meet the requirment for MC sim \n";
    }
    if (step%100 == 0) {
        PrintRate(step, true, true);
    }

} // for(int step=GetInitialStep(); step<GetFinalStep(); step++)
    std::clock_t end = std::clock();
    

    
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout<<"---- Simulation has ended ----\n";
    std::cout<<" The run took: "<<Nfunction::ConvertSecond2Time(elapsed_secs)<<"thread time and  ";
#ifdef _OPENMP
    double endwall_time = omp_get_wtime();
    endwall_time = endwall_time - startwall_time;
    std::cout<<Nfunction::ConvertSecond2Time(endwall_time)<<" wall time \n";
#endif

    m_pState->GetCurvatureCalculator()->Initialize();
    double Final_energy = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    double energy_leak = Final_energy - m_pState->GetEnergyCalculator()->GetEnergy();
    std::cout << std::fixed << std::setprecision(4);
    if(fabs(energy_leak) > 0.0001){
        
        std::cout<<"---> possible source of code error: energy leak... "<<energy_leak<<" with real energy of "<<Final_energy<<"  and stored energy of "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"\n";
    }

    
    
    double vol = 0;
    double g_c = 0;
    double t_a = 0;
    m_pState->GetVAHGlobalMeshProperties()->CalculateGlobalVariables(vol,t_a,g_c);
    vol -= m_pState->GetVAHGlobalMeshProperties()->GetTotalVolume();
    g_c -= m_pState->GetVAHGlobalMeshProperties()->GetTotalMeanCurvature();
    t_a -= m_pState->GetVAHGlobalMeshProperties()->GetTotalArea();

    if (m_pState->GetVAHGlobalMeshProperties()->VolumeIsActive())
    if(fabs(vol) > 0.0001 ){
        std::cout<<fabs(vol)<<" volume leak\n";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->GlobalCurvatureIsActive())
    if(fabs(g_c) > 0.0001)
    {
        std::cout<<fabs(g_c)<<" global curvature leak\n";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->AreaIsActive())
    if(fabs(t_a) > 0.0001)
    {
        std::cout<<fabs(t_a)<<" total area leak\n";
    }
        
    return true;
}
void  MC_Simulation::PrintRate(int step, bool clean, bool clear){
    
    double  vmove_rate =  100 * (m_pState->GetVertexPositionUpdate()->GetAcceptanceRate(clean));
    double  emove_rate =  100 * (m_pState->GetAlexanderMove()->GetAcceptanceRate(clean));
    double  imove_rate =  100 * (m_pState->GetInclusionPoseUpdate()->GetAcceptanceRate(clean));
    double  bmove_rate =  100 * (m_pState->GetDynamicBox()->GetAcceptanceRate(clean));
    double  vfmove_rate = 100 * (m_pState->GetVectorFieldsRotationUpdate()->GetAcceptanceRate(clean));

    std::cout<<"Step = "<<step<<"/"<<m_Final_Step<<std::flush;
    std::cout << std::fixed << std::setprecision(1);
    std::cout<<" Rates: "<<std::flush;
    std::cout<<" vertex move = "<<vmove_rate<<"%"<<std::flush;
    std::cout<<"; alexander move = "<<emove_rate<<"%"<<std::flush;
    if(m_pState->GetMesh()->GetInclusion().size() != 0){
        std::cout<<"; inclusion move = "<<imove_rate<<"%"<<std::flush;
    }
    if(m_pState->GetMesh()->GetNoVFPerVertex() != 0){
        std::cout<<"; vector fields move = "<<vfmove_rate<<"%"<<std::flush;
    }
    if(m_pState->GetDynamicBox()->GetDerivedDefaultReadName() != "No"){
        std::cout<<"; Box Move = "<<bmove_rate<<"%"<<std::flush;
    }
    if(clear){
        std::cout << '\r';
        std::cout << "\033[K";
    }
}
std::string MC_Simulation::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state + "\n Min_Max_Lenghts = "+Nfunction::D2S(m_MinLength2)+" "+Nfunction::D2S(m_MaxLength2);
    state = state + "\n MinfaceAngle = "+Nfunction::D2S(m_MinAngle);
    state = state + "\n Temprature = "+Nfunction::D2S(m_Beta)+" "+Nfunction::D2S(m_DBeta);
    state = state + "\n Box_Centering_F = "+Nfunction::D2S(m_CenteringFrequently);
    state = state + "\n Set_Steps = "+Nfunction::D2S(m_Initial_Step)+" "+Nfunction::D2S(m_Final_Step);
    
    return state;
}
bool MC_Simulation::CheckMesh(int step){
    
    if(m_CheckMeshFrequently == 0 || step%m_CheckMeshFrequently == 0){
        return true;
    }
    
    // this is only the  edges
    const std::vector<links *>& all_links = m_pState->GetMesh()->GetActiveL();
    for (std::vector<links *>::const_iterator it = all_links.begin() ; it != all_links.end(); ++it){

        vertex *p_v1 = (*it)->GetV1();
        vertex *p_v2 = (*it)->GetV2();
        double dist2 = p_v1->SquareDistanceFromAVertex(p_v2);
        if(dist2 < m_MinLength2 || dist2 > m_MaxLength2){
            
            return false;
        }
    }
   
    // all vertices should also have a distance larger then sqrt(m_MinLength2)


    
    return true;
}
