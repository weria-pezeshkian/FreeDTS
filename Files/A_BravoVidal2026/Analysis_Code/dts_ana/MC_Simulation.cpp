

#ifdef _OPENMP
# include <omp.h>
#endif
#ifdef MPI_DETECTED
# include <mpi.h>
#endif

#include <thread>
#include <ctime>
#include <iostream>
#include <complex>
#include <cmath>
#include <string.h>
#include <malloc.h>
#include "MC_Simulation.h"
#include "State.h"
#include "SimDef.h"
/*
 List of skipped function due to lack of clarity based on current state of the code
 They need to be finished before calling this a new version.
 
 1) void TimeSeriesLogInformation::WriteStateInfo(){
 
 */



using Complex = std::complex<double>;


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

#if MPI_DETECTED

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"Rank: "<<rank<<std::endl;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout<<"Size: "<<size<<std::endl;

#endif
    
//---> Voxelize the mesh for the first time; this should be done before any calculation

    
    
#if DEBUG_MODE == Enabled
    std::cout<<" system has been voxelaized  \n";
#endif

//----> checking if the mesh is good, within the bond of the simulation type. For here, it should be within
        //CheckMesh();
 
//--- ANA: No visualization!

   // time_t startTime;
   // time(&startTime);
#if DEBUG_MODE == Enabled
    std::cout<<" We have reached simulation run loop!  \n";
#endif

    
    
    std::vector<std::string> FramePath=m_pState->GetReadTrajTSI()->GetFilePaths();
    std::vector<int> FrameList=m_pState->GetReadTrajTSI()->GetFrameList();

    #ifdef MPI_DETECTED
    std::clock_t start = std::clock();
    if (rank==0){
    
    std::cout<<"------>   Simulation will be performed from "<<m_Initial_Step<<" to "<<m_Final_Step<<" steps\n";}
    #endif

    #ifndef MPI_DETECTED
    std::clock_t start = std::clock();
    std::cout<<"------>   Simulation will be performed from "<<m_Initial_Step<<" to "<<m_Final_Step<<" steps\n";
    #endif

    CreateMashBluePrint Create_BluePrint;
    //MeshBluePrint mesh_blueprint;


    if (m_pState->GetAnalysisVariables()->GetVisualizationActive()){
        m_pState->GetVisualization()->WriteAFrame(0);}

for (int step = m_Initial_Step; step < m_Final_Step; step++){

        

        std::string filename=FramePath[step];
        std::cout<<filename<<std::endl;

        /*std::ifstream input(m_pState->GetInputFile());
        ReadInclusionType(input,mesh);*/
        CreateMashBluePrint Create_BluePrint;
        MeshBluePrint mesh_blueprint = Create_BluePrint.MashBluePrintFromInput_Top(m_pState->GetInputFile(),filename);
        //std::ifstream input(m_pState->GetInputFile());
        //if (input.is_open()) {
        //    m_pState->ReadInclusionType(input);
        //    input.close();
        //    } else {
        //    std::cerr << "Error: Could not open input file for reading inclusion types." << std::endl;
        //}
        m_pState->GetMesh()->GenerateMesh(mesh_blueprint);
        
        

       
        
        //I need one specially thought for semi-flat membranees!!! A conditional is needed here
        if (m_pState->GetAnalysisVariables()->GetTopology()=="flat"){
            
            m_pState->GetMesh()->RemovePBC();
            m_pState->GetMesh()->CenterSemiFlatNOPBC();}
            //m_pState->GetMesh()->CenterNOPBC();}
            //m_pState->GetMesh()->CenterMesh();}
        else if (m_pState->GetAnalysisVariables()->GetTopology()=="spherical"){
            std::cout<<"Spherical Topology"<<std::endl;
             m_pState->GetMesh()->RemovePBC();
             m_pState->GetMesh()->CenterMeshSphericalNOPBC();
            }

        
        m_pState->GetCurvatureCalculator()->Initialize();
        
        m_pState->GetAnalysisCalculations()->Calculate();

        if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumActive()==true){
            if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumExtended()==0){
                m_pState->GetFluctationSpectrum()->CalculateSpectrum();}
            else if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumExtended()==1){
                m_pState->GetFluctationSpectrum()->CalculateSpectrumIndividual();
            }
            else if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumExtended()==2){
                m_pState->GetFluctationSpectrum()->CalculateSpectrumIndividualComplex();
            }
            else if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumExtended()==3){
                m_pState->GetFluctationSpectrum()->CalculateSpectrumSphericalComplex();
            }
        }

        if (m_pState->GetAnalysisVariables()->GetInclusionTrackingActive()==true){
            m_pState->GetInclusionTracking()->OutputInclusionInfo();
        }

        if (m_pState->GetAnalysisVariables()->GetInclusionClusterCalculationActive()==true){    
        
            m_pState->GetInclusionCluster()->CalculateClusterDistribution();}

        if (m_pState->GetAnalysisVariables()->GetRDFCalculationActive()==true){
            m_pState->GetRDF()->CalculateRDF();
        }
        
        m_pState->GetTimeSeriesDataOutput()->WriteTimeSeriesDataOutput(step);

        if (m_pState->GetAnalysisVariables()->GetVisualizationActive()){
            m_pState->GetVisualization()->WriteAFrame(step);}

        // Frames are read and the mesh is rebuilt from scratch every iteration,
        // which churns the heap; periodically release freed memory back to the OS
        // so RSS doesn't creep up over long trajectories.
        if ((step % 50) == 0){
            malloc_trim(0);
        }

} //End of simulation loop


    if (m_pState->GetAnalysisVariables()->GetFluctuationSpectrumActive()==true){
            m_pState->GetFluctationSpectrum()->CloseOutputStreams();
        }
// for(int step=GetInitialStep(); step<GetFinalStep(); step++)
     #ifdef MPI_DETECTED
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout<<"---- Simulation has ended ----\n";
    std::cout<<" The run took: "<<Nfunction::ConvertSecond2Time(elapsed_secs)<<"\n";}
    #endif

    #ifndef MPI_DETECTED
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    std::cout<<"---- Simulation has ended ----\n";
    std::cout<<" The run took: "<<Nfunction::ConvertSecond2Time(elapsed_secs)<<"\n";
    #endif 

    
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

bool MC_Simulation::ReadInclusionType(std::ifstream& input,MESH &mesh) {
    /*
     * @brief Reads and initializes inclusion types from an input file.
     *
     * This function reads the inclusion type definitions from an input file and initializes
     * the necessary data structures in the State object. It reads various properties for
     * each inclusion type and sets up the mesh to use these types.
     *
     * @param input Reference to the input file stream containing inclusion type definitions.
     * @return true if the inclusion types were successfully read and initialized, false otherwise.
     */
    
    std::string firstword, rest, str1, str2, TypeNames;
    int N, TypeID, NoType;
    double Kappa, KappaG, KappaP, KappaL, C0, C0P, C0N;

    // Store inclusion types in a vector
    std::vector<InclusionType> all_InclusionType;

    // Add a default inclusion type
    InclusionType emptyIncType;
    all_InclusionType.push_back(emptyIncType);

    // Read the header line
    input >> str1 >> NoType >> str2;
    getline(input, rest);
    getline(input, rest); // Discard the header line

    // Check if the header line indicates inclusion type definition
    if (str1 == "Define" || str1 == "define" || str1 == "DEFINE") {
        for (int i = 0; i < NoType; i++) {
            std::string inc_data;
            getline(input, inc_data);
            std::vector<std::string> inclusion_str = Nfunction::split(inc_data);
            N = Nfunction::String_to_Int(inclusion_str[0]);
            TypeNames = inclusion_str[1];
            Kappa = Nfunction::String_to_Double(inclusion_str[2]);
            KappaG = Nfunction::String_to_Double(inclusion_str[3]);
            KappaP = Nfunction::String_to_Double(inclusion_str[4]);
            KappaL = Nfunction::String_to_Double(inclusion_str[5]);
            C0 = Nfunction::String_to_Double(inclusion_str[6]);
            C0P = Nfunction::String_to_Double(inclusion_str[7]);
            C0N = Nfunction::String_to_Double(inclusion_str[8]);

            // Parse edge data if available
            double lam = 0, ekg = 0, ekn = 0, ecn = 0;
            if (inclusion_str.size() >= 13) {
                lam = Nfunction::String_to_Double(inclusion_str[9]);
                ekg = Nfunction::String_to_Double(inclusion_str[10]);
                ekn = Nfunction::String_to_Double(inclusion_str[11]);
                ecn = Nfunction::String_to_Double(inclusion_str[12]);
            }

            // Create inclusion type and add to vector
            InclusionType incType(TypeNames, i + 1, N, Kappa/2, KappaG, KappaP/2, KappaL/2, C0, C0P, C0N, lam, ekg, ekn, ecn);
            all_InclusionType.push_back(incType);
        }
    }

    // Set inclusion types in the mesh object
    mesh.m_InclusionType = all_InclusionType;

    // Set pointers to inclusion types
    mesh.m_pInclusionType.clear();
    for (size_t i = 0; i < mesh.m_InclusionType.size(); ++i) {
        mesh.m_pInclusionType.push_back(&mesh.m_InclusionType[i]);
    }

    return true;
}
