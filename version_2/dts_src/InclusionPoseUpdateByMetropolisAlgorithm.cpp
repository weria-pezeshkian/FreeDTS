

#include <stdio.h>
#include "State.h"
#include "InclusionPoseUpdateByMetropolisAlgorithm.h"

InclusionPoseUpdateByMetropolisAlgorithm::InclusionPoseUpdateByMetropolisAlgorithm (State *pState)  :
                m_pState(pState),
                m_pInclusion(pState->GetMesh()->GetInclusion()),
                m_Beta(pState->GetSimulation()->GetBeta()),
                m_DBeta(pState->GetSimulation()->GetDBeta())
{
                m_NumberOfMovePerStep_Angle  = 1;
                m_NumberOfMovePerStep_Kawasaki = 1;
}
InclusionPoseUpdateByMetropolisAlgorithm::InclusionPoseUpdateByMetropolisAlgorithm (State *pState, double rate_kawa, double rate_angle)  :
                    m_pState(pState),
                    m_pInclusion(pState->GetMesh()->GetInclusion()),
                    m_Beta(pState->GetSimulation()->GetBeta()),
                    m_DBeta(pState->GetSimulation()->GetDBeta())
{
                    m_NumberOfMovePerStep_Angle = rate_kawa;
                    m_NumberOfMovePerStep_Kawasaki = rate_angle;
}
InclusionPoseUpdateByMetropolisAlgorithm::~InclusionPoseUpdateByMetropolisAlgorithm() {
    
}
bool InclusionPoseUpdateByMetropolisAlgorithm::Initialize() {
    return true;
}
bool InclusionPoseUpdateByMetropolisAlgorithm::EvolveOneStep(int step){
    
    int no_incs = m_pInclusion.size();
    int no_steps_angle = no_incs*m_NumberOfMovePerStep_Angle;
    int no_steps_kawa = no_incs* m_NumberOfMovePerStep_Kawasaki;
    
  for (int i = 0; i< no_steps_kawa;i++) {
    
      int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(no_incs);
      inclusion *p_inc = m_pInclusion[r_inc_id];
      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      
      std::vector <links*> linklist = p_inc->Getvertex()->GetVLinkList();
      int m = m_pState->GetRandomNumberGenerator()->IntRNG(linklist.size());
      links *d_link = linklist[m];
      
      if(KawasakiMove(p_inc,d_link, thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }
    for (int i = 0; i< no_steps_angle;i++) {
      
        int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(no_incs);
        inclusion *p_inc = m_pInclusion[r_inc_id];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        
        if(RotationMove(p_inc, m_DR*dx, thermal)){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
      }
    
    return true;
}
bool InclusionPoseUpdateByMetropolisAlgorithm::KawasakiMove(inclusion* p_inc, links * d_links, double temp) {
    /**
     * @brief Perform a Kawasaki move on the inclusion using the Metropolis algorithm.
     *
     * This function attempts to move an inclusion from its current vertex to a target vertex
     * while preserving its local pose. It uses the Metropolis algorithm to decide whether to accept the move
     * based on the change in energy and a given temperature.
     *
     * @param p_inc Pointer to the inclusion to be moved.
     * @param d_links Pointer to the link connecting the source and target vertices.
     * @param temp The temperature parameter for the Metropolis acceptance criterion.
     * @return true if the move is accepted, false otherwise.
     */
    
    double old_energy = 0.0;
    double new_energy = 0.0;

    vertex * hver = p_inc->Getvertex();      // Current vertex of the inclusion
    vertex * tver = d_links->GetV2();        // Target vertex for the inclusion
    bool has_inclusion = tver->VertexOwnInclusion();   // Check if target vertex has an inclusion
    
    
    // If both vertices contain the same type of inclusion, no move is needed
   if (has_inclusion && p_inc->GetInclusionType() == tver->GetInclusion()->GetInclusionType()) {
        return false;
    }

    double en_1 = hver->GetEnergy();
    double en_2 = tver->GetEnergy();
    old_energy = en_1 + en_2;
    
    // add the interaction energies and make a copy to the links associated with them
    std::vector<links*> Affected_links = GetEdgesWithInteractionChange(d_links);
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        old_energy += 2*(*it)->GetIntEnergy();
    }

//--> performing the moves
    if(has_inclusion) {
        // give the host vertex the traget vertex inclusion (only is allowed if target vertex has one)
        hver->UpdateInclusion((tver->GetInclusion()));
        (tver->GetInclusion())->Updatevertex(hver);
    }
    else { // if the second vertex does not have an inclusion, the first one will have none now
         hver->UpdateOwnInclusion(false);
    }
    // move the specific inc
    p_inc->Updatevertex(tver);
    tver->UpdateInclusion(p_inc);
    tver->UpdateOwnInclusion(true);

    // Calculate new energy
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(hver);
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(tver);

    //-- interaction energy should be calculated here
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
    }
    

    //--> elatsic energy
    double diff_energy = new_energy - old_energy;

    double U = m_Beta * diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > temp ) {
        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        return true;
    }
    else {
//---->  Move rejected, revert changes
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
        }
        //---> revert the vertices energy
        hver->UpdateEnergy(en_1);
        tver->UpdateEnergy(en_2);

        if(has_inclusion) {
            tver->UpdateInclusion(hver->GetInclusion());
            tver->GetInclusion()->Updatevertex(tver);
        }
        else { // if the second vertex did not have an inclusion, it should return to that
            tver->UpdateOwnInclusion(false);
            hver->UpdateOwnInclusion(true);
        }
        p_inc->Updatevertex(hver);
        hver->UpdateInclusion(p_inc);

        return false;
     }
    return true;

}
//bool InclusionPoseUpdateByMetropolisAlgorithm::RotationMove(inclusion *p_inc, double dx, double dy, double temp) {
bool InclusionPoseUpdateByMetropolisAlgorithm::RotationMove(inclusion *p_inc, double dx, double temp) {

    /**
     * @brief Perform a rotational move on the inclusion using the Metropolis algorithm.
     *
     * This function attempts to rotate an inclusion while preserving its local pose and
     * calculates the resulting energy change to decide whether to accept the move based on the Metropolis criterion.
     *
     * @param p_inc Pointer to the inclusion to be rotated.
     * @param dx Change in the x-direction of the local direction vector.
     * @param dy Change in the y-direction of the local direction vector.
     * @param temp The temperature parameter for the Metropolis acceptance criterion.
     * @return true if the move is accepted, false otherwise.
     */
    
    double old_energy = 0.0;
    double new_energy = 0.0;
    
    vertex * ver = p_inc->Getvertex();
    std::vector<links*> Affected_links = ver->GetVLinkList();
    if(ver->GetVertexType() == 1){
        Affected_links.push_back(ver->GetPrecedingEdgeLink());
    }
    double en_1 = ver->GetEnergy();
    old_energy += en_1;
    
    // Copy interaction energies and get old energy
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
    }
    
//-- perfrom the move
     Vec3D local_d = p_inc->GetLDirection();
     Vec3D Dr(local_d(1), -local_d(0), 0);
     Dr = local_d + (Dr * dx);
     Dr.normalize();
     p_inc->UpdateLocalDirection(Dr);
     
    
    
    
    new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(ver);

    //-- interaction energy should be calculated here
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
    }
    //--> elatsic energy
    double diff_energy = new_energy - old_energy;
    
    double U = m_Beta * diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > temp ) {
        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        return true;
    }
    else {
//---->  Move rejected, revert changes
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
        }
        //---> revert the vertices energy
        ver->UpdateEnergy(en_1);
        p_inc->UpdateLocalDirection(local_d);  // orginal orinetation

        return false;
     }

    return true;
}
std::vector<links*> InclusionPoseUpdateByMetropolisAlgorithm::GetEdgesWithInteractionChange(links* p_link){
    /**
     * @brief Get the list of links that  potentially changing interaction energies due to inclsuion swich.
     *
     * This function identifies and returns a list of links that might have their interaction
     * energies affected. It includes neighboring links of the vertices connected by `p_link` and ensures no duplicates.
     * @param p_link Pointer to the link whose interaction change is being considered.
     * @return std::vector<links*> A list of links with potentially changing interaction energies.
     */
    /*
               \ /    \ /
                v1 === v2
               / \     / \
     */
    std::vector<links *> Affected_links; // return list
    
    vertex * ver_1 = p_link->GetV1();
    vertex * ver_2 = p_link->GetV2();
    
    // Collect neighboring links for both vertices
    std::vector<links *> neighbour_link1 = ver_1->GetVLinkList();     // links that their interaction energies may change
    std::vector<links *> neighbour_link2 = ver_2->GetVLinkList();     // links that their interaction energies may change
    
    // Insert all neighboring links into the set
     Affected_links.insert(Affected_links.end(), neighbour_link1.begin(), neighbour_link1.end());   // a copy of the links
     Affected_links.insert(Affected_links.end(), neighbour_link2.begin(), neighbour_link2.end());   // a copy of the links

    // Insert preceding edge links if vertices are of type 1
    if(ver_1->GetVertexType() == 1)
        Affected_links.push_back(ver_1->GetPrecedingEdgeLink());
    if(ver_2->GetVertexType() == 1)
        Affected_links.push_back(ver_2->GetPrecedingEdgeLink());

    // Erase the original link and its mirror link if necessary
    // note, it may look silly, but when we have edge links, is hard to know what do remove ....
    
    Affected_links.erase(std::remove(Affected_links.begin(), Affected_links.end(), p_link), Affected_links.end());
    if(p_link->GetMirrorFlag()){
        Affected_links.erase(std::remove(Affected_links.begin(), Affected_links.end(), p_link->GetMirrorLink()), Affected_links.end());
    }
    // now added back the original link (not sure)
    Affected_links.push_back(p_link);
    
    return Affected_links;

}
std::string InclusionPoseUpdateByMetropolisAlgorithm::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Kawasaki) +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Angle);
    return state;
}









