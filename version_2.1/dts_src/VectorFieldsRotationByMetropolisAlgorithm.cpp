

#include <stdio.h>
#include "State.h"
#include "VectorFieldsRotationByMetropolisAlgorithm.h"

VectorFieldsRotationByMetropolisAlgorithm::VectorFieldsRotationByMetropolisAlgorithm (State *pState)  :
                m_pState(pState),
                m_pActiveV(pState->GetMesh()->GetActiveV()),
                m_Beta(pState->GetSimulation()->GetBeta()),
                m_DBeta(pState->GetSimulation()->GetDBeta()),
                m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex())
{
    m_DR = 0.1;
    m_NumberOfMovePerStep  = 1;
}
VectorFieldsRotationByMetropolisAlgorithm::VectorFieldsRotationByMetropolisAlgorithm (State *pState, double rate)  :
                    m_pState(pState),
                    m_pActiveV(pState->GetMesh()->GetActiveV()),
                    m_Beta(pState->GetSimulation()->GetBeta()),
                    m_DBeta(pState->GetSimulation()->GetDBeta()),
                    m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex())
{
                    m_DR = 0.1;
                    m_NumberOfMovePerStep = rate;
}
VectorFieldsRotationByMetropolisAlgorithm::VectorFieldsRotationByMetropolisAlgorithm (State *pState, double rate, double dr)  :
                    m_pState(pState),
                    m_pActiveV(pState->GetMesh()->GetActiveV()),
                    m_Beta(pState->GetSimulation()->GetBeta()),
                    m_DBeta(pState->GetSimulation()->GetDBeta()),
                    m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex())
{
                    m_DR = dr;
                    m_NumberOfMovePerStep = rate;
}
VectorFieldsRotationByMetropolisAlgorithm::~VectorFieldsRotationByMetropolisAlgorithm() {
    
}
bool VectorFieldsRotationByMetropolisAlgorithm::Initialize() {
    return true;
}
bool VectorFieldsRotationByMetropolisAlgorithm::EvolveOneStep(int step){
    
    if(m_No_VectorFields_Per_V == 0){
        return false;
    }
    int no_of_V = m_pActiveV.size();
    int no_steps = no_of_V*m_NumberOfMovePerStep;
    

    for (int layer = 0; layer< m_No_VectorFields_Per_V; layer++) {
        for (int i = 0; i< no_steps; i++) {
            int r_v_id = m_pState->GetRandomNumberGenerator()->IntRNG(no_of_V);
            double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            double dx = 1 - 2 * (m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        
            if(RotationMove(layer, m_pActiveV[r_v_id], m_DR*dx, thermal)){
                m_AcceptedMoves++;
            }
            m_NumberOfAttemptedMoves++;
      }
    }
    
    return true;
}
bool VectorFieldsRotationByMetropolisAlgorithm::RotationMove(int layer, vertex *p_vertex, double dx, double temp){
    
    /**
     * @brief Perform a rotational move a vector field using the Metropolis algorithm.
     *
     * This function attempts to rotate an inclusion while preserving its local pose and
     * calculates the resulting energy change to decide whether to accept the move based on the Metropolis criterion.
     */
    
    double old_energy = 0.0;
    double new_energy = 0.0;
    VectorField *p_vf = p_vertex->GetVectorField(layer);

    
    std::vector<links*> Affected_links = p_vertex->GetVLinkList();
    if(p_vertex->GetVertexType() == 1){
        Affected_links.push_back(p_vertex->GetPrecedingEdgeLink());
    }
     double en_1 = p_vf->GetMembraneBindingEnergy(); //
     old_energy += en_1;
    
    // Copy interaction energies and get old energy
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
    
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy(layer);
    }
    
//-- perfrom the move
    Vec3D local_d = p_vf->GetLDirection();
    Vec3D Dr(local_d(1), -local_d(0), 0);
    Dr = local_d + (Dr * dx);
    Dr.normalize();
    p_vf->UpdateLocalDirection(Dr);
    
    double en_2 = (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(p_vf,p_vertex);
    new_energy += en_2;

    //-- interaction energy should be calculated here
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        
          new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(layer, *it);
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
        p_vf->UpdateMembraneBindingEnergy(en_1);
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_VFInteractionEnergy();
        }
        //---> revert the vertices energy
        p_vf->UpdateLocalDirection(local_d);  // orginal orinetation

        return false;
     }

    return true;
}
std::string VectorFieldsRotationByMetropolisAlgorithm::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep) +" "+ Nfunction::D2S(m_DR);
    return state;
}









