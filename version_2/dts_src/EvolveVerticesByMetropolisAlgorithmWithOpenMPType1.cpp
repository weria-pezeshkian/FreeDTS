

#include <stdio.h>
#include <thread>   // For std::this_thread::sleep_for
#include <chrono>   // For std::chrono::milliseconds
#ifdef _OPENMP
# include <omp.h>
#endif
#include "EvolveVerticesByMetropolisAlgorithmWithOpenMPType1.h"
#include "State.h"

EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState)
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
EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState, double rate_surf, double rate_edge, double dr)
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
EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::~EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(){

}
void EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::Initialize(){
    m_pBox = m_pState->GetMesh()->GetBox();
    m_Total_ThreadsNo = m_pState->GetThreads_Number();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    return;
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveOneStep(int step){
 
#ifdef _OPENMP
    // OpenMP is available, continue with the program
    
    int no_surf_v = m_pSurfV.size();
    int no_edge_v = m_pEdgeV.size();
    int no_steps_edge = no_edge_v*m_NumberOfMovePerStep_Edge;
    int no_steps_surf = no_surf_v*m_NumberOfMovePerStep_Surf;




    //omp_set_num_threads(m_Total_ThreadsNo);  // Set the number of threads
    double diff_energy = 0.0;
    int total_AcceptedMoves = 0;

    // Use OpenMP to parallelize the outer loop with dynamic scheduling and reduction
    #pragma omp parallel reduction(+:total_AcceptedMoves, diff_energy)
    {
        int local_AcceptedMoves = 0;  // Thread-local counters
        double local_diff_energy = 0.0;
        vertex* pvertex;

        // Parallel for with dynamic scheduling
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < no_steps_surf; i++) {

            // Implement a lock with backoff to avoid busy-waiting
            while (true) {
                int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v);  // Thread-local RNG
                pvertex = m_pSurfV[r_vid];

                // Try to lock the vertex and its neighbors
                if (pvertex->CheckLockVertex()) {
                    if (!pvertex->CheckLockNeighbourVertex()) {
                        pvertex->UnlockVertex();  // Unlock if neighbors can't be locked
                    } else {
                        break;  // Successfully locked vertex and neighbors, break the loop
                    }
                }
            }
            /*int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v); // Thread-local RNG
            pvertex = m_pSurfV[r_vid];
        while (true) {
            if (pvertex->CheckLockVertex()) {
                if (!pvertex->CheckLockNeighbourVertex()) {
                    pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
                } else {
                    break;  // Successfully locked vertex and neighbors, break the loop
                }
            }
        std::this_thread::sleep_for(std::chrono::nanoseconds(Backoff_Factor));
        }*/

            // Skip if the vertex is in the frozen group
            if (m_FreezGroupName == pvertex->GetGroupName()) {
                pvertex->UnlockVertex();
                pvertex->UnlockNeighbourVertex();
                continue;
            }

            // Generate random move vectors
            double dx = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double dy = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double dz = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

            // Check if the move is within the boundary
            if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx, dy, dz, pvertex)) {
                pvertex->UnlockVertex();
                pvertex->UnlockNeighbourVertex();
                continue;
            }

            // Generate random thermal value
            double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double tem_en = 0.0;

            // Evolve the vertex and update local counters if successful
            if (EvolveOneVertex(step, pvertex, m_DR * dx, m_DR * dy, m_DR * dz, thermal, tem_en)) {
                local_AcceptedMoves++;  // Increment local accepted moves
                local_diff_energy += tem_en;  // Accumulate local energy change
            }

            // Unlock the vertex and its neighbors
            pvertex->UnlockVertex();
            pvertex->UnlockNeighbourVertex();
        }

        // Update the shared variables after the loop
        total_AcceptedMoves += local_AcceptedMoves;  // Reduction for accepted moves
        diff_energy += local_diff_energy;  // Reduction for energy difference
    }

    // Update the global variables once, after the parallel section
    m_AcceptedMoves += total_AcceptedMoves;  // Safely update accepted moves
    m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);  // Add energy
    m_NumberOfAttemptedMoves += no_steps_surf;  // Update the number of attempted moves

omp_set_num_threads(m_Total_ThreadsNo); // Set the number of threads
 diff_energy = 0.0;

// Use OpenMP to parallelize the outer loop with dynamic scheduling
#pragma omp parallel for schedule(dynamic) num_threads(m_Total_ThreadsNo)
for (int i = 0; i < no_steps_edge; i++) {
    // Each thread has its own local counters
    int local_AcceptedMoves = 0;
    double local_diff_energy = 0.0;

    vertex *pvertex;

    // Try to lock a random vertex and its neighbors
    /*int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v); // Thread-local RNG
    	pvertex = m_pEdgeV[r_vid];
    while (true) {
        if (pvertex->CheckLockVertex()) {
            if (!pvertex->CheckLockNeighbourVertex()) {
                pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
            } else {
                break;  // Successfully locked vertex and neighbors, break the loop
            }
        }
    std::this_thread::sleep_for(std::chrono::nanoseconds(Backoff_Factor));
    }*/
    
    
    while (true) {
        int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v);
        pvertex = m_pEdgeV[r_vid];

        if (pvertex->CheckLockVertex()) {
            if (!pvertex->CheckLockNeighbourVertex()) {
                pvertex->UnlockVertex();  // Unlock the vertex if neighbors can't be locked
            } else {
                break;  // Successfully locked vertex and neighbors, break the loop
            }
        }
    }

    // Skip if vertex is in the frozen group
    if (m_FreezGroupName == pvertex->GetGroupName()) {
        pvertex->UnlockVertex();
        pvertex->UnlockNeighbourVertex();
        continue;
    }

    // Generate random move vectors
    double dx = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double dy = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double dz = 1 - 2 * m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    // Check if the move is within the boundary
    if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx, dy, dz, pvertex)) {
        pvertex->UnlockVertex();
        pvertex->UnlockNeighbourVertex();
        continue;
    }

    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    double tem_en = 0.0;

    // Evolve the vertex and update local counters if successful
    if (EvolveOneVertex(step, pvertex, m_DR * dx, m_DR * dy, m_DR * dz, thermal, tem_en)) {
        local_AcceptedMoves++;
        local_diff_energy += tem_en;
    }

    // Unlock the vertex and its neighbors
    pvertex->UnlockVertex();
    pvertex->UnlockNeighbourVertex();

    // Combine local results into shared variables
    #pragma omp atomic
    m_AcceptedMoves += local_AcceptedMoves;

    // Atomic update to shared diff_energy variable
    #pragma omp atomic
    diff_energy += local_diff_energy;
}

// Add total energy from the local diff_energy calculated by each thread
m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
m_NumberOfAttemptedMoves += no_steps_surf;



    return true;
#else
    std::cerr << "---> Error: OpenMP is not available, but the program requires it. Please recompile with OpenMP support.\n";
    exit(EXIT_FAILURE); // Exit with a non-zero value to indicate failure
    
    return false;
#endif
}
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp, double &changed_en){
    
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
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2){
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
        double dist2 = (*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
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
                    if((*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it) < mindist2)
                        return false;
                }
                
            }
        } ///   for(int s=-1;s<2;s++){
    return true;
}
// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::CheckFacesAfterAVertexMove(vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(m_MinAngle) || !link->CheckFaceAngleWithNextEdgeFace(m_MinAngle)) {
            return false;
        }
    }
    return true;
}

//========
// this function can be deleted any time; it is for test cases only
double  EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::SystemEnergy()
{
    /*
    MESH* m_pMESH = m_pState->m_pMesh;
    std::vector<vertex *> ActiveV = m_pMESH->m_pActiveV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = mLink.begin() ; it != mLink.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    
    for (std::vector<vertex *>::iterator it = ActiveV.begin() ; it != ActiveV.end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
    en=m_pEnergyCalculator->TotalEnergy(ActiveV,mLink);
    //en=en+ m_pEnergyCalculator->TotalEnergy(EdgeV,pEdgeL);
   
    
    return en;
     */
    return 0;
}
std::string EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::CurrentState(){
    
    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Surf) +" "+Nfunction::D2S(m_NumberOfMovePerStep_Edge);
    state = state +" "+ Nfunction::D2S(m_DR);
    return state;
}
std::vector<links*> EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::GetEdgesWithInteractionChange(vertex* p_vertex){
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



