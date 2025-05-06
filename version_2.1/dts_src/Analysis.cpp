

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include <ctime>
#include <iostream>
#include <string.h>
#include "Analysis.h"
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
Analysis::Analysis(State *pState, std::string path) : m_pState(pState), m_Path(path) {

    
}
Analysis::~Analysis(){
    
}
void Analysis::Initialize(){
    
    return;
}
bool Analysis::do_Simulation(){

//---> Voxelize the mesh for the first time; this should be done before any calculation
    
    std::string InputFileName = m_pState->GetInputFileName();
    std::string gname = m_pState->GetRunTag();
    if (m_Path.back() != '/') {  // Ensure trailing slash
        m_Path += "/";
    }
    std::ofstream edgeFrame;
    edgeFrame.open("edge.xvg");
for (int step = m_Initial_Step; step <= m_Final_Step; step++){
        
    std::string TopologyFile = m_Path + gname + Nfunction::Int_to_String(step)+"."+TSIExt;
    if(!Nfunction::FileExist(TopologyFile)){
        std::cout<<"--> warning: "<<TopologyFile<<" does not exist, we skip it \n";
        continue;
    }
    CreateMashBluePrint Create_BluePrint;
    MeshBluePrint mesh_blueprint = Create_BluePrint.MashBluePrintFromInput_Top(InputFileName, TopologyFile);
    m_pState->GetMesh()->Clear_Mesh();
    m_pState->GetMesh()->GenerateMesh(mesh_blueprint);
    m_pState->GetCurvatureCalculator()->Initialize();
    CenterMesh();
    std::vector<Vec3D > ActiveForce;
    const std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    Vec3D tot_F(0,0,0);
    for (std::vector<vertex *>::const_iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it){
        Vec3D Fr = m_pState->GetForceonVerticesfromInclusions()->Inclusion_Force((*it));
        tot_F = tot_F + Fr;
        ActiveForce.push_back(Fr);
    }
    std::cout<<tot_F<<"\n";
    int period = m_pState->GetVisualization()->GetPeriod();
    m_pState->GetVisualization()->ClearVector();
    m_pState->GetVisualization()->AddVector("active", ActiveForce);
    m_pState->GetVisualization()->WriteAFrame(step*period);
    //m_pState->GetVoxelization()->SetBox(m_pMesh->GetBox());
    //m_pState->GetVoxelization()->Voxelize(m_pState->GetMesh()->GetActiveV());
    
    
    const std::vector<vertex *>& edge_vertex = m_pState->GetMesh()->GetEdgeV();
    edgeFrame<<"frame "<<step<<"  "<<edge_vertex.size()<<"\n";
    if (!edge_vertex.empty()) {
        vertex * oneedgeV = edge_vertex[0];
        edgeFrame << oneedgeV->GetPos() << "\n";
        vertex * nextV = oneedgeV->GetEdgeLink()->GetV2();
        while (true) {
            if (nextV != oneedgeV) {
                edgeFrame << nextV->GetPos() << "\n";
                nextV = nextV->GetEdgeLink()->GetV2();
            } else {
                break;
            }
        }
    }

    
    std::cout<<step<<" step \n";
    

    
    
} // for(int step=GetInitialStep(); step<GetFinalStep(); step++)

    return true;
}
std::string Analysis::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();

    
    return state;
}
bool Analysis::CheckMesh(int step){
    
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
bool Analysis::CenterMesh(){
    
    
    Vec3D COG(0,0,0);
    Vec3D *pBox = m_pState->GetMesh()->GetBox();
    const std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    double tot_V = double(all_vertex.size());
    for (std::vector<vertex *>::const_iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it){
        COG = COG + (*it)->GetPos()*(1.0/tot_V);
    }
    for (std::vector<vertex *>::const_iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it){
        Vec3D newPos = (*it)->GetPos()- COG + (*pBox)*(0.5);
        (*it)->UpdateVXPos(newPos(0));
        (*it)->UpdateVYPos(newPos(1));
        (*it)->UpdateVZPos(newPos(2));

    }
    
    return true;
}
