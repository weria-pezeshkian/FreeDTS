

#include <stdio.h>
#include <thread>   // For std::this_thread::sleep_for
#include <chrono>   // For std::chrono::milliseconds
#ifdef _OPENMP
# include <omp.h>
#endif
#include "EvolveVerticesByMetropolisAlgorithmWithOpenMPType2.h"
#include "State.h"
#include "./Registry/FactoryVertexPositionIntegrator.h"

EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::EvolveVerticesByMetropolisAlgorithmWithOpenMPType2(State *pState)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetDBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()){
          
          m_NumberOfMovePerStep_Surf = 1.0;
          m_NumberOfMovePerStep_Edge = 1.0;
          m_DR = 0.05;
          
      }
EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::EvolveVerticesByMetropolisAlgorithmWithOpenMPType2(State *pState, double rate_surf, double rate_edge, double dr)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetDBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()) {
          
          m_NumberOfMovePerStep_Surf = rate_surf;
          m_NumberOfMovePerStep_Edge = rate_edge;
          m_DR = dr;
      }
EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::~EvolveVerticesByMetropolisAlgorithmWithOpenMPType2(){

}
void EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::Initialize(){
    m_pBox = m_pState->GetMesh()->GetBox();
    m_Total_ThreadsNo = m_pState->GetThreads_Number();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    return;
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::EvolveOneStep(int step) {
#ifdef _OPENMP
    const int no_surf_v = m_pSurfV.size();
    const int no_edge_v = m_pEdgeV.size();
    const int no_steps_edge = no_edge_v * m_NumberOfMovePerStep_Edge;
    const int no_steps_surf = no_surf_v * m_NumberOfMovePerStep_Surf;
    
    // Process surface vertices
    double diff_energy = ProcessVertices(step, m_pSurfV, no_steps_surf);
    
    // Process edge vertices
    diff_energy += ProcessVertices(step, m_pEdgeV, no_steps_edge);
    
    return true;
#else
    std::cerr << "---> Error: OpenMP is not available, but the program requires it. "
              << "Please recompile with OpenMP support.\n";
    exit(EXIT_FAILURE);
    return false;
#endif
}

double EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::ProcessVertices(
    int step, const std::vector<vertex*>& vertices, int num_steps) {
    double total_diff_energy = 0.0;
#ifdef _OPENMP
    const int no_vertices = vertices.size();
    int total_accepted_moves = 0;
    
    // Use reduction for all cases - OpenMP optimizes this well
    #pragma omp parallel reduction(+:total_accepted_moves, total_diff_energy)
    {
        int local_accepted = 0;
        double local_diff_energy = 0.0;
        
        // Each thread gets its own RNG for thread safety
        auto* RNG = m_pState->GetRandomNumberGenerator();
        
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < num_steps; ++i) {
            ProcessSingleVertexAttempt(step, vertices, no_vertices, RNG,
                                      local_accepted, local_diff_energy);
        }
        
        total_accepted_moves += local_accepted;
        total_diff_energy += local_diff_energy;
    }
    
    // Update global state
    m_AcceptedMoves += total_accepted_moves;
    m_pState->GetEnergyCalculator()->AddToTotalEnergy(total_diff_energy);
    m_NumberOfAttemptedMoves += num_steps;
    
#endif
    return total_diff_energy;

}

void EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::ProcessSingleVertexAttempt(
    int step, const std::vector<vertex*>& vertices,
    int no_vertices, RNG* rng,
    int& local_accepted, double& local_diff_energy) {
#ifdef _OPENMP

    
    vertex* pvertex = nullptr;
    bool lock_acquired = false;
    
    // Lock vertex with backoff strategy
    const int max_retries = 100;
    for (int retry = 0; retry < max_retries; ++retry) {
        const int r_vid = rng->IntRNG(no_vertices);
        pvertex = vertices[r_vid];
        
        if (TryLockVertexWithNeighbors(pvertex)) {
            lock_acquired = true;
            break;  // Successfully locked
        }
        
        // Backoff: yield after 10 attempts to reduce contention
        if (retry >= 10) {
            std::this_thread::yield();
        }
    }
    
    // Failed to acquire lock after retries - skip this iteration
    if (!lock_acquired) {
        return;
    }
    
    // Check frozen group
    if (m_FreezGroupName == pvertex->GetGroupName()) {
        UnlockVertexWithNeighbors(pvertex);
        return;
    }
    
    // Generate random move vector
    const double dx = 1.0 - 2.0 * rng->UniformRNG(1.0);
    const double dy = 1.0 - 2.0 * rng->UniformRNG(1.0);
    const double dz = 1.0 - 2.0 * rng->UniformRNG(1.0);
    
    // Check boundary
    if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(
            m_DR * dx, m_DR * dy, m_DR * dz, pvertex)) {
        UnlockVertexWithNeighbors(pvertex);
        return;
    }
    
    // Attempt move
    const double thermal = rng->UniformRNG(1.0);
    double tem_en = 0.0;
    
    if (EvolveOneVertex(step, pvertex, m_DR * dx, m_DR * dy,
                       m_DR * dz, thermal, tem_en)) {
        ++local_accepted;
        local_diff_energy += tem_en;
    }
    
    UnlockVertexWithNeighbors(pvertex);
#endif

}

bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::TryLockVertexWithNeighbors(
    vertex* pvertex) {
#ifdef _OPENMP
    if (!pvertex->CheckLockVertex()) {
        return false;
    }
    
    if (!pvertex->CheckLockNeighbourVertex()) {
        pvertex->UnlockVertex();
        return false;
    }
#endif

    return true;
}

void EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::UnlockVertexWithNeighbors(
    vertex* pvertex) {
#ifdef _OPENMP
    pvertex->UnlockVertex();
    pvertex->UnlockNeighbourVertex();
#endif
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp, double &changed_en){
    
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    if(!VertexMoveIsFine(pvertex,dx,dy,dz,m_MinLength2,m_MaxLength2))  // this function could get a booling varaible to say, it crossed the voxel
        return 0;

    //--- obtain vertices energy terms and make copies
    old_energy = pvertex->GetEnergy();
    old_energy += pvertex->GetBindingEnergy();
    pvertex->ConstantMesh_Copy();
    pvertex->Copy_VFsBindingEnergy();  // vector field
    const std::vector<vertex *>& vNeighbourV = pvertex->GetVNeighbourVertex();  
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        (*it)->ConstantMesh_Copy();
        old_energy += (*it)->GetEnergy();
        old_energy += (*it)->GetBindingEnergy();
        (*it)->Copy_VFsBindingEnergy();


    }
    std::vector<triangle *> N_triangles = pvertex->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
        (*it)->ConstantMesh_Copy();
    }
    const std::vector<links *>& v_NLinks = pvertex->GetVLinkList();
    for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
        
        (*it)->ConstantMesh_Copy();
        (*it)->GetNeighborLink1()->ConstantMesh_Copy();
    }
    //-- we need this to make sure all the links connected to this v is updated
    if(pvertex->GetVertexType() == 1){
        pvertex->GetPrecedingEdgeLink()->ConstantMesh_Copy();
    }
    // find the links in which there interaction energy changes
    std::vector<links*> Affected_links = GetEdgesWithInteractionChange(pvertex);
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }
    // --- obtaining global variables that can change by the move. Note, this is not the total volume, only the one that can change.
     double old_Tvolume = 0;
     double old_Tarea = 0;
     double old_Tcurvature = 0;
     double new_Tvolume = 0;
     double new_Tarea = 0;
     double new_Tcurvature = 0;
//--->
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, old_Tvolume, old_Tarea, old_Tcurvature);
    }
    //---> for now, only active nematic force: ForceonVerticesfromInclusions
    Vec3D Dx(dx,dy,dz);
    double dE_force_from_inc  = m_pState->GetForceonVerticesfromInclusions()->Energy_of_Force(pvertex, Dx);
    double dE_force_from_vector_fields  = m_pState->GetForceonVerticesfromVectorFields()->Energy_of_Force(pvertex, Dx);
    double dE_force_on_vertex  = m_pState->GetForceonVertices()->Energy_of_Force(pvertex, Dx);

//----> Move the vertex;
        pvertex->PositionPlus(dx,dy,dz);
    //--- update triangles normal
    for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
        (*it)->UpdateNormal_Area(m_pBox);
    }
    //  check new faces angles, if bad, reverse the trinagles
    if(!CheckFacesAfterAVertexMove(pvertex)){
        
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        pvertex->PositionPlus(-dx,-dy,-dz);
        return false;
    }
//---->
    //--> calculate edge shape operator;
    for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
        
       // (*it)->UpdateNormal();
        //  (*it)->GetNeighborLink1()->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
        (*it)->GetNeighborLink1()->UpdateShapeOperator(m_pBox);
    }
    //-- we need this to make sure all the links connected to this v is updated
    if(pvertex->GetVertexType() == 1){
        pvertex->GetPrecedingEdgeLink()->UpdateEdgeVector(m_pBox);
    }
    // --> calculate vertex shape operator
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(pvertex);
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(*it);
    }
    //---> calculate new energies
    new_energy = (m_pState->GetEnergyCalculator())->SingleVertexEnergy(pvertex);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(pvertex);

    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(*it);
        new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(*it);

    }
    //-- interaction energy should be calculated here

    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        
        if(pvertex->GetNumberOfVF() != 0 ){
            for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                
                new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
        }
    }
    //---> get energy for ApplyConstraintBetweenGroups
    double dE_Cgroup = m_pState->GetApplyConstraintBetweenGroups()->CalculateEnergyChange(pvertex, Dx);

//---> new global variables
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
    
    //--> only elatsic energy
    double diff_energy = new_energy - old_energy;
            changed_en = diff_energy;
    //std::cout<<diff_energy<<" dif en \n";
    //--> sum of all the energies
    double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + dE_force_from_inc + dE_force_from_vector_fields + dE_volume + dE_t_area + dE_g_curv ;
    double U = m_Beta * tot_diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > temp ) {
        // move is accepted
     // (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);        
        //---> if vertex is out of the voxel, update its voxel
        if(!pvertex->CheckVoxel()){
            pvertex->UpdateVoxelAfterAVertexMove();
        }
        //---> ApplyConstraintBetweenGroups
        m_pState->GetApplyConstraintBetweenGroups()->AcceptMove();
        
        //---> global variables
        if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
           
            m_pState->GetVAHGlobalMeshProperties()->Add2Volume(new_Tvolume - old_Tvolume);
            m_pState->GetVAHGlobalMeshProperties()->Add2TotalArea(new_Tarea - old_Tarea);
            m_pState->GetVAHGlobalMeshProperties()->Add2GlobalCurvature(new_Tcurvature - old_Tcurvature);
        }
        return true;
    }
    else {
//---> reverse the changes that has been made to the system
        //---> reverse the triangles
        changed_en = 0;
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        //---> reverse the links
        for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
            
            (*it)->ReverseConstantMesh_Copy();
            (*it)->GetNeighborLink1()->ReverseConstantMesh_Copy();
        }
        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
            (*it)->Reverse_VFInteractionEnergy();

        }
        //-- we need this to make sure all the links connected to this v is updated
        if(pvertex->GetVertexType() == 1){
            pvertex->GetPrecedingEdgeLink()->ReverseConstantMesh_Copy();
        }
        //---> reverse the vertices
        pvertex->ReverseConstantMesh_Copy();
        pvertex->Reverse_VFsBindingEnergy();

        for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
            (*it)->Reverse_VFsBindingEnergy();
        }

        return false;
     }
    return true;
}
//---> this does not check the angle of the faces. Because this should be done after the move:
//waste of calculation if we do ith before the move. Unless, we store the values. That also not good because move could get rejected.
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2){
//--->  vertex new position if get accepted
    
        double new_x = pvertex->GetXPos() + dx;
        double new_y = pvertex->GetYPos() + dy;
        double new_z = pvertex->GetZPos() + dz;
        //-- if the adding crosses the box
        if (new_x >= (*m_pBox)(0)) {
            new_x = new_x - (*m_pBox)(0);
        } else if (new_x < 0) {
            new_x = new_x + (*m_pBox)(0);
        }
        if (new_y >= (*m_pBox)(1)) {
            new_y = new_y - (*m_pBox)(1);
        } else if (new_y < 0) {
            new_y = new_y + (*m_pBox)(1);
        }
        if (new_z >= (*m_pBox)(2)) {
            new_z = new_z - (*m_pBox)(2);
        } else if (new_z < 0) {
            new_z = new_z + (*m_pBox)(2);
        }
  
//--->  let check the distances with the nighbours
    std::vector <vertex *> npvertex = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = npvertex.begin() ; it != npvertex.end(); ++it){
        double dist2 = pvertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
            if(dist2 < mindist2 || dist2 > maxdist2)
            return false;
    }
//---> now check it within the voxel cells
//---> lets get the object voxel, note: the voxel could be different from the associated one as it is moving
        //-- obtain the object (vertex) new cell id, with respect to the current cell
        int NoX = pvertex->GetVoxel()->GetXNoVoxel();
        int NoY = pvertex->GetVoxel()->GetYNoVoxel();
        int NoZ = pvertex->GetVoxel()->GetZNoVoxel();
    
        int new_nx = int(new_x/pvertex->GetVoxel()->GetXSideVoxel((*m_pBox)(0)));
        int new_ny = int(new_y/pvertex->GetVoxel()->GetYSideVoxel((*m_pBox)(1)));
        int new_nz = int(new_z/pvertex->GetVoxel()->GetZSideVoxel((*m_pBox)(2)));
    
        int old_nx = pvertex->GetVoxel()->GetXIndex();
        int old_ny = pvertex->GetVoxel()->GetYIndex();
        int old_nz = pvertex->GetVoxel()->GetZIndex();
    
        int i = Voxel<int>::Convert2LocalVoxelIndex(new_nx, old_nx, NoX);
        int j = Voxel<int>::Convert2LocalVoxelIndex(new_ny, old_ny, NoY);
        int k = Voxel<int>::Convert2LocalVoxelIndex(new_nz, old_nz, NoZ);

    
        //-- check if it has moved too far
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            std::cout << " ---> warning: the object might moved more than one voxel " << std::endl;
            
#if DEBUG_MODE == Enabled
std::cout << i<<" "<<j<<"  "<<k<<"  local "<< std::endl;
std::cout << int(new_x/pvertex->GetVoxel()->GetXSideVoxel((*m_pBox)(0)))<<" "<<int(new_y/pvertex->GetVoxel()->GetYSideVoxel((*m_pBox)(1)))<<"  "<<int(new_z/pvertex->GetVoxel()->GetZSideVoxel((*m_pBox)(2)))<<"  new voxel "<< std::endl;
std::cout << pvertex->GetVoxel()->GetXIndex()<<" "<<pvertex->GetVoxel()->GetYIndex()<<"  "<<pvertex->GetVoxel()->GetZIndex()<<" old voxel"<< std::endl;
#endif
            return false;
        }
  
        Voxel<vertex>* new_pvox = pvertex->GetVoxel()->GetANeighbourCell(i, j, k);
  
        for(int n=-1;n<2;n++)
        for(int m=-1;m<2;m++)
        for(int s=-1;s<2;s++){
            std::vector <vertex *> CV = new_pvox->GetANeighbourCell(n, m, s)->GetContentObjects();
            for (std::vector<vertex *>::iterator it = CV.begin() ; it != CV.end(); ++it){
                if(*it != pvertex){
                    if(pvertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it) < mindist2)
                        return false;
                }
                
            }
        } ///   for(int s=-1;s<2;s++){
    ///
    ///
    /*  */
    return true;
}
// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::CheckFacesAfterAVertexMove(vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(m_MinAngle) || !link->CheckFaceAngleWithNextEdgeFace(m_MinAngle)) {
            return false;
        }
    }
    return true;
}


std::string EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::CurrentState(){
    
    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Surf) +" "+Nfunction::D2S(m_NumberOfMovePerStep_Edge);
    state = state +" "+ Nfunction::D2S(m_DR);
    return state;
}
std::vector<links*> EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::GetEdgesWithInteractionChange(vertex* p_vertex){
#if DEVELOPMENT_MODE == Enabled
    std::cout<<" DEVELOPMENT_MODE ID 665656: This function should be made much better\n ";
    // for example, precheck to select only links that both has includioon ...
#endif
    
    std::vector<links*> edge_with_interaction_change;
    std::vector<links *> all_temlinks;
    
    std::vector<vertex *> neighbor_vertices = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = neighbor_vertices.begin() ; it != neighbor_vertices.end(); ++it)
    {
            if((*it)->VertexOwnInclusion() || p_vertex->GetNumberOfVF() != 0 ) {  // due to vector fields
                std::vector<links *> ltem = (*it)->GetVLinkList();
                all_temlinks.insert(all_temlinks.end(), ltem.begin(), ltem.end());
                // note, this is even need it if the p_vertex vertex is not an edge vertex
                if((*it)->m_VertexType == 1){
                    all_temlinks.push_back((*it)->m_pPrecedingEdgeLink);
                }
            }
    }

    
    //-- now we remove the repeated links
    for (std::vector<links *>::iterator it = all_temlinks.begin() ; it != all_temlinks.end(); ++it)
    {
        bool Should_be_added = true;
        for (std::vector<links *>::iterator it2 = edge_with_interaction_change.begin() ; it2 != edge_with_interaction_change.end(); ++it2)
        {
            if( *it2 == *it){
                Should_be_added = false;
            }
            else if((*it2)->GetMirrorFlag() && (*it2)->GetMirrorLink() == *it){
                Should_be_added = false;
            }
        }
        if(Should_be_added == true)
            edge_with_interaction_change.push_back((*it));

    }
    
    return edge_with_interaction_change;
}
namespace
{
#ifdef _OPENMP

    AbstractVertexPositionIntegrator* Create_VerPosMC_OpenMPType2(
        std::istream& input,
        State* state)
    {
        double rate_surf = 0.0;
        double rate_edge = 0.0;
        double dr = 0.0;

        if (!(input >> rate_surf >> rate_edge >> dr)) {
            return nullptr;
        }

        std::string rest;
        std::getline(input, rest);

        return new EvolveVerticesByMetropolisAlgorithmWithOpenMPType2(
            state,
            rate_surf,
            rate_edge,
            dr
        );
    }

#else

    // Stub: OpenMP not available
    AbstractVertexPositionIntegrator* Create_VerPosMC_OpenMPType2(
        std::istream&,
        State*)
    {
        std::cerr
            << "Error: '"
            << EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::GetDefaultReadName()
            << "' requires OpenMP, but this build does not support it.\n"
            << "' you can just change to the single thread algorithm .\n";

        return nullptr;
    }

#endif

    const bool registered = []()
    {
        FactoryVertexPositionIntegrator::Instance().Register(
            EvolveVerticesByMetropolisAlgorithmWithOpenMPType2::GetDefaultReadName(),
            &Create_VerPosMC_OpenMPType2
        );
        return true;
    }();
}

