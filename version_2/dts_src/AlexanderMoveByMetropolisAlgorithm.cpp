

#include <stdio.h>
#include "AlexanderMoveByMetropolisAlgorithm.h"
#include "State.h"
AlexanderMoveByMetropolisAlgorithm::AlexanderMoveByMetropolisAlgorithm(State* pState) :
            m_pState(pState),
            m_pSurfL(pState->GetMesh()->GetRightL()),
            m_Beta(pState->GetSimulation()->GetBeta()),
            m_DBeta(pState->GetSimulation()->GetDBeta()),
            m_MinLength2(pState->GetSimulation()->GetMinL2()),
            m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
            m_MinAngle(pState->GetSimulation()->GetMinAngle()) {
                
                m_NumberOfMovePerStep = 1;
}
AlexanderMoveByMetropolisAlgorithm::AlexanderMoveByMetropolisAlgorithm(State* pState, double rate) :
            m_pState(pState),
            m_pSurfL(pState->GetMesh()->GetRightL()),
            m_Beta(pState->GetSimulation()->GetBeta()),
            m_DBeta(pState->GetSimulation()->GetDBeta()),
            m_MinLength2(pState->GetSimulation()->GetMinL2()),
            m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
            m_MinAngle(pState->GetSimulation()->GetMinAngle()){
                
                m_NumberOfMovePerStep = rate;

}
AlexanderMoveByMetropolisAlgorithm::~AlexanderMoveByMetropolisAlgorithm(){
    
}
bool AlexanderMoveByMetropolisAlgorithm::Initialize() {
    
    m_pBox = m_pState->GetMesh()->GetBox();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::EvolveOneStep(int step){
 
    int no_edges = m_pSurfL.size();
    int no_steps = no_edges*m_NumberOfMovePerStep;
    
  for (int i = 0; i< no_steps;i++) {
    
      int r_lid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edges);
      links *p_link = m_pSurfL[r_lid];

      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      if(FlipOneEdge(step, p_link,thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }
    
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::FlipOneEdge(int step, links *p_edge, double temp){
    /**
     * @brief Attempts to flip the specified edge and evaluates the energy change.
     *
     * This function attempts to flip the given edge in the mesh. It checks if the edge can be flipped,
     * calculates the energy before and after the flip, and decides whether to accept the flip based on
     * the Metropolis criterion. If the flip is not accepted, it reverts the changes.
     *
     * @param step The current simulation step.
     * @param p_edge A pointer to the edge to be flipped.
     * @param temp The current temperature of the system for the Metropolis criterion.
     * @return true if the edge flip was accepted, false otherwise.
     */
    
    double old_energy = 0.0;
    double new_energy = 0.0;

    // Check if the edge can be flipped; if not, return false
    if (!EdgeCanBeFliped(p_edge)) {
        return false;
    }

    // Obtain references to the vertices connected by the edge and its mirror
    vertex *v1 = p_edge->GetV1();
    vertex *v2 = p_edge->GetV2();
    vertex *v3 = p_edge->GetV3();
    vertex *v4 = p_edge->GetMirrorLink()->GetV3();

    old_energy = v1->GetEnergy();
    old_energy += v2->GetEnergy();
    old_energy += v3->GetEnergy();
    old_energy += v4->GetEnergy();
    
    old_energy += v1->GetBindingEnergy();
    old_energy += v2->GetBindingEnergy();
    old_energy += v3->GetBindingEnergy();
    old_energy += v4->GetBindingEnergy();


//-- get the energy for interaction
    std::vector<links*> Affected_links = GetEdgesWithInteractionChange(p_edge);
    
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();

    }
    
    // Obtain and sum the initial global variables that might change
    double old_Tvolume = 0.0, old_Tarea = 0.0, old_Tcurvature = 0.0;
    double new_Tvolume = 0.0, new_Tarea = 0.0, new_Tcurvature = 0.0;
//--->
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateALinkTrianglesContributionToGlobalVariables(p_edge, old_Tvolume, old_Tarea, old_Tcurvature);
    }
//---> make a copy of the edge

    links *ml  = p_edge->GetMirrorLink();
    links *l1  = p_edge->GetNeighborLink1();
    links *l2  = p_edge->GetNeighborLink2();
    links *l3  = ml->GetNeighborLink1();
    links *l4  = ml->GetNeighborLink2();

    triangle *t1 = p_edge->GetTriangle();
    triangle *t2 = ml->GetTriangle();

//--- make a copy of the triangles
    t1->ConstantMesh_Copy();
    t2->ConstantMesh_Copy();

//-- fliping the link
    if(!p_edge->Flip(p_edge)) {
        return false;
    }
    
//--- update the triangle normal
    t1->UpdateNormal_Area(m_pBox);
    t2->UpdateNormal_Area(m_pBox);
    
// Check if the resulting faces are valid
    if (!CheckFacesAfterFlip(p_edge)) {
        p_edge->Reverse_Flip(p_edge);
        t1->ReverseConstantMesh_Copy();
        t2->ReverseConstantMesh_Copy();
        return false;
    }
// we proceed with the changes and calculate energies
    // make copy. since we reverse the flip by reverse flip function, we only need to copy with constant mesh, meaning no list will be changed
    p_edge->ConstantMesh_Copy();
    l1->ConstantMesh_Copy();
    l2->ConstantMesh_Copy();
    l3->ConstantMesh_Copy();
    l4->ConstantMesh_Copy();
    
    v1->ConstantMesh_Copy();
    v2->ConstantMesh_Copy();
    v3->ConstantMesh_Copy();
    v4->ConstantMesh_Copy();
    v1->Copy_VFsBindingEnergy();
    v2->Copy_VFsBindingEnergy();
    v3->Copy_VFsBindingEnergy();
    v4->Copy_VFsBindingEnergy();

    

//-- Update geometry
    //--- update the edges normal
    p_edge->UpdateShapeOperator(m_pBox);
    l1->UpdateShapeOperator(m_pBox);
    l2->UpdateShapeOperator(m_pBox);
    l3->UpdateShapeOperator(m_pBox);
    l4->UpdateShapeOperator(m_pBox);

    // --> calculate vertex shape operator
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(v1);
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(v2);
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(v3);
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(v4);

//--- obtain vertices energy terms
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(v1);
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(v2);
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(v3);
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(v4);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v1);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v2);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v3);
    new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v4);

//-- interaction energy should be calculated here
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
          new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
    }
//---> new global variables
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateALinkTrianglesContributionToGlobalVariables(p_edge, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
    
    double diff_energy = new_energy - old_energy;
    //--> sum of all the energies
    double tot_diff_energy = diff_energy +  dE_volume + dE_t_area + dE_g_curv ;

    double U = m_Beta * tot_diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > temp ) {
        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        
        //---> global variables
        if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
             m_pState->GetVolumeCoupling()->UpdateArea_Volume(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
            m_pState->GetTotalAreaCoupling()->UpdateTotalArea(old_Tarea, new_Tarea);
            m_pState->GetGlobalCurvature()->UpdateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
        }
        return true;
    }
    else {
//---> reverse the changes that has been made to the system
        p_edge->Reverse_Flip(p_edge);

        //---> reverse the triangles
        t1->ReverseConstantMesh_Copy();
        t2->ReverseConstantMesh_Copy();
        
        p_edge->ReverseConstantMesh_Copy();
        l1->ReverseConstantMesh_Copy();
        l2->ReverseConstantMesh_Copy();
        l3->ReverseConstantMesh_Copy();
        l4->ReverseConstantMesh_Copy();

        //---> reverse the vertices
        v1->ReverseConstantMesh_Copy();
        v2->ReverseConstantMesh_Copy();
        v3->ReverseConstantMesh_Copy();
        v4->ReverseConstantMesh_Copy();
        v1->Reverse_VFsBindingEnergy();
        v2->Reverse_VFsBindingEnergy();
        v3->Reverse_VFsBindingEnergy();
        v4->Reverse_VFsBindingEnergy();

        
        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
            (*it)->Reverse_VFInteractionEnergy();
        }
        return false;

     }
    return true;
}
// this function can be deleted any time; it is for test cases only
double  AlexanderMoveByMetropolisAlgorithm::SystemEnergy()
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


// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool AlexanderMoveByMetropolisAlgorithm::CheckFacesAfterFlip(links* p_links) {
    // Checks the validity of faces after an edge flip operation
    // Parameters:
    //   p_links: Pointer to the edge to be checked
    // Returns:
    //   True if the faces remain valid after the flip, false otherwise
    
    links* p_mlinks = p_links->GetMirrorLink();
    const Vec3D& n1 = p_links->GetTriangle()->GetNormalVector();
    const Vec3D& n2 = p_mlinks->GetTriangle()->GetNormalVector();
    
    if( Vec3D::dot(n1,n2) < m_MinAngle)
        return false;
    
    if((p_links->GetNeighborLink1())->m_LinkType == 0) {
        
        const Vec3D& n11 = p_links->GetNeighborLink1()->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n1,n11) < m_MinAngle)
            return false;
    }
    
    if((p_links->GetNeighborLink2())->m_LinkType == 0) {
        
        const Vec3D& n12 = p_links->GetNeighborLink2()->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n1,n12) < m_MinAngle)
            return false;
    }
    if((p_mlinks->GetNeighborLink1())->m_LinkType == 0) {
        
        const Vec3D& n21 = p_mlinks->GetNeighborLink1()->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n2,n21) < m_MinAngle)
            return false;
    }
    if((p_mlinks->GetNeighborLink2())->m_LinkType == 0) {
        
        const Vec3D& n22 = p_mlinks->GetNeighborLink2()->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n2,n22) < m_MinAngle)
            return false;
    }
    
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::EdgeCanBeFliped(links *p_edge) {
    /**
     * @brief Checks if an edge can be flipped.
     *
     * It first checks if the vertices at the ends of the edge have
     * more than three neighbors. If not, the edge cannot be flipped. Then it ensures that
     * the vertex opposite to one end of the edge is not a neighbor of the vertex opposite
     * the other end of the edge. Additionally, it verifies that the distance between these
     * two vertices, considering periodic boundary conditions, is within the specified
     * allowable range.
     *
     * @param p_edge A pointer to the edge (of type 'links') that is being evaluated.
     * @return true if the edge can be flipped, false otherwise.
     *
     * @details
     * The function performs the following steps:
     * 1. Obtains the neighbor counts for the vertices at the ends of the edge.
     * 2. Checks if these neighbor counts are greater than three.
     * 3. Fetches the neighbors of one of the vertices.
     * 4. Ensures the opposite vertex is not already a neighbor.
     * 5. Computes the distance between the two vertices with periodic boundary conditions.
     * 6. Verifies that this distance is within the acceptable range.
     *
     * The function is optimized for performance by reducing repeated calls to getters
     * and using references where applicable.
     */
    
// Obtain references to the vertices connected by the edge
    vertex* v1 = p_edge->GetV1();
    vertex* v2 = p_edge->GetV2();
    vertex* v3 = p_edge->GetV3();
    vertex* v4 = p_edge->GetMirrorLink()->GetV3();

    // Get the sizes of the neighbor vertices for v1 and v2
    int v1NeighborSize = v1->GetVNeighbourVertex().size();
    int v2NeighborSize = v2->GetVNeighbourVertex().size();

    // If either vertex has 3 or fewer neighbors, the edge cannot be flipped
    if (v1NeighborSize <= 3 || v2NeighborSize <= 3)
        return false;

    // Get the neighbor vertices of v3 as a reference
    const std::vector<vertex *> &v3Neighbors = v3->GetVNeighbourVertex();

    // Check if v4 is already a neighbor of v3
    for (std::vector<vertex *>::const_iterator it = v3Neighbors.begin(); it != v3Neighbors.end(); ++it) {
        if (*it == v4)
            return false;
    }

    // Calculate the distance between v3 and v4 considering periodic boundary conditions
    double dx = v4->GetVXPos() - v3->GetVXPos();
    if (fabs(dx) > (*m_pBox)(0) / 2.0) {
        if (dx < 0)
            dx += (*m_pBox)(0);
        else
            dx -= (*m_pBox)(0);
    }

    double dy = v4->GetVYPos() - v3->GetVYPos();
    if (fabs(dy) > (*m_pBox)(1) / 2.0) {
        if (dy < 0)
            dy += (*m_pBox)(1);
        else
            dy -= (*m_pBox)(1);
    }

    double dz = v4->GetVZPos() - v3->GetVZPos();
    if (fabs(dz) > (*m_pBox)(2) / 2.0) {
        if (dz < 0)
            dz += (*m_pBox)(2);
        else
            dz -= (*m_pBox)(2);
    }

    double dist_2 = dx * dx + dy * dy + dz * dz;

    // Check if the distance squared is within the allowed range
    if (dist_2 > m_MaxLength2 || dist_2 < m_MinLength2)
        return false;

    return true;
}
std::vector<links*> AlexanderMoveByMetropolisAlgorithm::GetEdgesWithInteractionChange(links *p_edge){
    /**
     * @brief Retrieves a list of edges whose interactions may change due to the flipping of a given edge.
     *
     * This function determines which edges' interactions might be affected if the specified edge is flipped.
     * It checks the vertices of the edge and collects all edges connected to these vertices that have certain
     * properties, ensuring no duplicates are added to the result.
     *
     * @param p_edge A pointer to the edge (of type 'links') that is being evaluated.
     * @return A vector of pointers to edges (of type 'links') that may have interaction changes.
     *
     * @details
     * The function performs the following steps:
     * 1. Fetches the vertices connected by the edge and their opposing vertices.
     * 2. Collects all edges connected to these vertices, taking into account their inclusion status and type.
     * 3. Ensures that edges are not duplicated in the result list.
     */
    
// Fetch vertices connected by the edge and their opposing vertices
    vertex* v1 = p_edge->GetV1();
    vertex* v2 = p_edge->GetV2();
    vertex* v3 = p_edge->GetV3();
    vertex* v4 = p_edge->GetMirrorLink()->GetV3();
    
// Vector to store edges that may have interaction changes
    std::vector<links*> edge_with_interaction_change;
// Temporary list for gathering edges from the vertices
    std::vector <links*> temlinklist;
    
// Gather edges from v1 if it owns an inclusion
        if(v1->VertexOwnInclusion() || v1->GetNumberOfVF() != 0  )
        {
            
            std::vector<links *> ltem = v1->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            
            if(v1->m_VertexType==1)
            temlinklist.push_back(v1->m_pPrecedingEdgeLink);
        }
        if(v2->VertexOwnInclusion() || v2->GetNumberOfVF() != 0  )
        {
            
            std::vector<links *> ltem=v2->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(v2->m_VertexType==1)
            temlinklist.push_back(v2->m_pPrecedingEdgeLink);
        }
        if(v3->VertexOwnInclusion() || v3->GetNumberOfVF() != 0  )
        {
            
            std::vector<links *> ltem = v3->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(v3->m_VertexType==1)
            temlinklist.push_back(v3->m_pPrecedingEdgeLink);
        }
        if(v4->VertexOwnInclusion() || v4->GetNumberOfVF() != 0  )
        {
            std::vector<links *> ltem=v4->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(v4->m_VertexType==1)
            temlinklist.push_back(v4->m_pPrecedingEdgeLink);
        }

// Remove duplicates and mirror edges from the list
    for (std::vector<links *>::iterator it = temlinklist.begin() ; it != temlinklist.end(); ++it)
    {
        bool Should_be_added = true;
        for (std::vector<links *>::iterator it2 = edge_with_interaction_change.begin() ; it2 != edge_with_interaction_change.end(); ++it2)
        {
            if((*it2) == (*it)){
                
                Should_be_added = false;
            }
            else if((*it2)->GetMirrorFlag() && ((*it2)->GetMirrorLink())==(*it)){
                
                Should_be_added = false;
            }
        }
        if(Should_be_added)
            edge_with_interaction_change.push_back((*it));
    }
    
    return edge_with_interaction_change;
}
std::string AlexanderMoveByMetropolisAlgorithm::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + " "+ Nfunction::D2S(m_NumberOfMovePerStep);
    return state;
}
