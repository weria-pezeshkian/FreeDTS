#include "OpenEdgeEvolutionWithConstantVertex.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */

OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex(int period, double rate, State *pState ) :
                                                    m_pState(pState),
                                                    m_pEdgeL(pState->GetMesh()->GetEdgeL()),
                                                    m_pGhostL(pState->GetMesh()->GetGhostL()),
                                                    m_pGhostT(pState->GetMesh()->GetGhostT()),
                                                    m_pActiveT(pState->GetMesh()->GetActiveT()),
                                                    m_pRightL(pState->GetMesh()->GetRightL()),
                                                    m_pLeftL(pState->GetMesh()->GetLeftL()),
                                                    m_pActiveL(pState->GetMesh()->GetActiveL()),
                                                    m_pSurfV(pState->GetMesh()->GetSurfV()),
                                                    m_pEdgeV(pState->GetMesh()->GetEdgeV()),
                                                    m_Period(period),
                                                    m_NumberOfMovePerStep(rate),
                                                    m_Beta(pState->GetSimulation()->GetBeta()),
                                                    m_DBeta(pState->GetSimulation()->GetDBeta()),
                                                    m_MinLength2(pState->GetSimulation()->GetMinL2()),
                                                    m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
                                                    m_MinAngle(pState->GetSimulation()->GetMinAngle()),
                                                    m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex())
{
    m_EdgeSize = m_pEdgeL.size();
}
OpenEdgeEvolutionWithConstantVertex::~OpenEdgeEvolutionWithConstantVertex(){
    
}
void OpenEdgeEvolutionWithConstantVertex::Initialize(){
    
    m_pBox = m_pState->GetMesh()->GetBox();
    m_pMesh = m_pState->GetMesh();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;
    
    m_EdgeSize = m_pEdgeL.size();

    return;
}
bool OpenEdgeEvolutionWithConstantVertex::Move(int step) {

    if(step % m_Period != 0 || m_pEdgeL.size() == 0 )
        return false;
   
    int N = m_pEdgeL.size();
    N = int(m_NumberOfMovePerStep*double(N));
    for (int i = 0; i< N ;i++) {

        if(m_pEdgeL.size() == 0 )
            break;
        if(MCAttemptedToAddALink()){
            m_AcceptedMoves++;
        }
        if(MCAttemptedToRemoveALink()){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
        m_NumberOfAttemptedMoves++;
    }
    m_EdgeSize = m_pEdgeL.size();

    return true;
}
bool OpenEdgeEvolutionWithConstantVertex::MCAttemptedToRemoveALink(){
    
    
    if( m_pEdgeL.size() == 0 )
        return false;
    
    int n = m_pState->GetRandomNumberGenerator()->IntRNG(m_pEdgeL.size());
    links *plink = m_pEdgeL[n];
    
    vertex *v1 = plink->GetV1();
    vertex *v2 = plink->GetV2();
    vertex *v3 = plink->GetV3();

    if(v1->GetVLinkList().size() < 2 || v2->GetVLinkList().size() < 2 || v3->GetVertexType() == 1)
        return false;
        
    double eold = 0;
    double enew = 0;
    
        // calculate the old energy
        eold = v1->GetEnergy();
        eold += v2->GetEnergy();
        eold += v3->GetEnergy();
        eold += v1->GetBindingEnergy();
        eold += v2->GetBindingEnergy();
        eold += v3->GetBindingEnergy();

        //
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
                eold += 2 * (*it)->GetIntEnergy();
                eold += 2 * (*it)->GetVFIntEnergy();
            }
        
        std::vector <links *> nvl2 = v2->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
                eold += 2 * (*it)->GetIntEnergy();
                eold += 2 * (*it)->GetVFIntEnergy();
            }
        
        std::vector <links *> nvl3 = v3->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
                eold += 2 * (*it)->GetIntEnergy();
                eold += 2 * (*it)->GetVFIntEnergy();
            }
            
            //--- this link does not exist in the nb of any vertices
            eold += 2 * (v1->m_pPrecedingEdgeLink)->GetIntEnergy();
            eold += 2 * (v1->m_pPrecedingEdgeLink)->GetVFIntEnergy();

            // because these two links were counted two times
            eold -= 2 * (plink->GetNeighborLink1())->GetIntEnergy();
            eold -= 2 * (plink->GetNeighborLink2())->GetIntEnergy();
            eold -= 2 * (plink->GetNeighborLink1())->GetVFIntEnergy();
            eold -= 2 * (plink->GetNeighborLink2())->GetVFIntEnergy();

        }
        //== this is a new piece to use global inclusion direction as a constant variable during the move: note global direction should be projected on the new local plane
        if(v1->VertexOwnInclusion()){
            if(!v1->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }
        if(v2->VertexOwnInclusion()){
            if(!v2->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }
        if(v3->VertexOwnInclusion()){
            if(!v3->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }
        // check if the g_direction for vector field can be found
        if(!v1->UpdateVFGlobalDirectionFromLocalDirection()){
            return false;
        }
        if(!v2->UpdateVFGlobalDirectionFromLocalDirection()){
            return false;
        }
        if(!v3->UpdateVFGlobalDirectionFromLocalDirection()){
            return false;
        }
        // we kill a link and update the geomotry
        KillALink(plink);

        // updating the local inc direction from the global one
        bool inc_direction = true;
        bool vf_direction = true;
            if(v1->VertexOwnInclusion()){
                inc_direction = v1->GetInclusion()->UpdateLocalDirectionFromGlobal();
            }
            if(v2->VertexOwnInclusion()){
                inc_direction = v2->GetInclusion()->UpdateLocalDirectionFromGlobal();
            }
            if(v3->VertexOwnInclusion()){
                inc_direction = v3->GetInclusion()->UpdateLocalDirectionFromGlobal();
            }
        // update local from global for vector fields
        vf_direction = v1->UpdateVFLocalDirectionFromGlobalDirection();
        vf_direction = v2->UpdateVFLocalDirectionFromGlobalDirection();
        vf_direction = v3->UpdateVFLocalDirectionFromGlobalDirection();
    
    // if local could not be updated from global
        if(!vf_direction || !inc_direction){
            v1->ReverseVFLocalDirection();
            v2->ReverseVFLocalDirection();
            v3->ReverseVFLocalDirection();
            if(v1->VertexOwnInclusion()){
                v1->GetInclusion()->Reverse_Direction();
            }
            if(v2->VertexOwnInclusion()){
                v2->GetInclusion()->Reverse_Direction();
            }
            if(v3->VertexOwnInclusion()){
                v3->GetInclusion()->Reverse_Direction();
            }
            CreateALink(v1);
            return false;
        }

    
        // new energy
            enew = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
            enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
            enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
            enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v1);
            enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v2);
            enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v3);


        
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
                enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
            
        std::vector <links *> nvl2 = v2->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
                enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
            
        std::vector <links *> nvl3 = v3->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
                enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
        }
        }
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->GetPrecedingEdgeLink());
        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, v1->GetPrecedingEdgeLink());
        }
        
        double diff_energy = enew - eold;
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    
        //if(double(NS)/double(NL+1)*(exp(-m_Beta*DE)>thermal ))
    if(exp( -m_Beta * diff_energy + m_DBeta) > thermal) {

            m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
            return true;
        }
        else{
            // reject the move
                CreateALink(v1);
            v1->ReverseVFLocalDirection();
            v2->ReverseVFLocalDirection();
            v3->ReverseVFLocalDirection();
            if(v1->VertexOwnInclusion()){
                v1->GetInclusion()->Reverse_Direction();
            }
            if(v2->VertexOwnInclusion()){
                v2->GetInclusion()->Reverse_Direction();
            }
            if(v3->VertexOwnInclusion()){
                v3->GetInclusion()->Reverse_Direction();
            }
            
            {
            double e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
            e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
            e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);

            e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v1);
            e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v2);
            e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v3);

            
                std::vector <links *> nvl1 = v1->GetVLinkList();
                    for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
                         m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                            (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                        }
                    }
                    
                std::vector <links *> nvl2 = v2->GetVLinkList();
                    for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
                        m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                            (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                        }
                    }
                    
                std::vector <links *> nvl3 = v3->GetVLinkList();
                for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
                    m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                    for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                        (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                    }
                }
                m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->GetPrecedingEdgeLink());
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, v1->GetPrecedingEdgeLink());
                }
            }
            return false;

           // std::cout<<"killed rejected \n";

        }
        // int energy;
        //double eint = TwoInclusionsInteractionEnergy(l0);
    
    return true; // MCAttemptedToRemoveALink(){
} // MCAttemptedToRemoveALink(){

bool OpenEdgeEvolutionWithConstantVertex::MCAttemptedToAddALink(){
    
    // an atempt to create a link                       //       v3
                                                        //      /   \
                                                      //     v1-----v2
    
    
    bool ClosingAHole = false; // a flag to check if v2-1 link exist, so this means the chosen hole will be closed
    
    if(m_pEdgeL.size() < 3 && m_pEdgeL.size() != 0){
        std::cout<<" error--> 123 should not happen edge number is "<<m_pEdgeL.size()<<"\n";
        return false;
    }

    // select an edge vertex
    int n = m_pState->GetRandomNumberGenerator()->IntRNG(m_pEdgeV.size());
    vertex *v1 = m_pEdgeV[n];
    
    if( !Linkisvalid(v1) ){
        return false;
    }
    
    double eold = 0;
    double enew = 0;
    
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    
    // if this is true, means that the hole will be closed
    if(v2->m_pEdgeLink->GetV2() == v1){
        ClosingAHole = true;
    }
    
    
    
//----   for constant global direction type of moves
    if(v1->VertexOwnInclusion()){
        if(!(v1->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    if(v2->VertexOwnInclusion()){
        if(!(v2->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    if(v3->VertexOwnInclusion()){
        if(!(v3->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    // check if the g_direction for vector field can be found
    if(!v1->UpdateVFGlobalDirectionFromLocalDirection()){
        return false;
    }
    if(!v2->UpdateVFGlobalDirectionFromLocalDirection()){
        return false;
    }
    if(!v3->UpdateVFGlobalDirectionFromLocalDirection()){
        return false;
    }
    eold += v1->GetEnergy();
    eold += v2->GetEnergy();
    eold += v3->GetEnergy();
    
    eold += v1->GetBindingEnergy();
    eold += v2->GetBindingEnergy();
    eold += v3->GetBindingEnergy();

    //=== inclusion interaction energy; it is multipled by two because only half is assinged to each edge
    std::vector <links *> nvl1 = v1->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
        eold += 2 * (*it)->GetIntEnergy();
        eold += 2 * (*it)->GetVFIntEnergy();

    }

    std::vector <links *> nvl2 = v2->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
        eold += 2 * (*it)->GetIntEnergy();
        eold += 2 * (*it)->GetVFIntEnergy();
    }

    std::vector <links *> nvl3 = v3->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
        eold += 2 * (*it)->GetIntEnergy();
        eold += 2 * (*it)->GetVFIntEnergy();
    }

    eold += 2 * (v1->m_pPrecedingEdgeLink)->GetIntEnergy();
    eold += 2 * (v1->m_pPrecedingEdgeLink)->GetVFIntEnergy();

    // create a link (this also updates the gemotry)
    // Only one of these two will be valid
    links *newlink;
    triangle *T;
    
    if(ClosingAHole){
        T = CloseATriangleHole(v1);
    }
    else{
        newlink = CreateALink(v1);

    }
  //  return true;
    //----   for constant global direction type of moves
    // updating the local inc direction from the global one
    bool inc_direction = true;
    bool vf_direction = true;
        if(v1->VertexOwnInclusion()){
            inc_direction = v1->GetInclusion()->UpdateLocalDirectionFromGlobal();
        }
        if(v2->VertexOwnInclusion()){
            inc_direction = v2->GetInclusion()->UpdateLocalDirectionFromGlobal();
        }
        if(v3->VertexOwnInclusion()){
            inc_direction = v3->GetInclusion()->UpdateLocalDirectionFromGlobal();
        }
    // update local from global for vector fields
    vf_direction = v1->UpdateVFLocalDirectionFromGlobalDirection();
    vf_direction = v2->UpdateVFLocalDirectionFromGlobalDirection();
    vf_direction = v3->UpdateVFLocalDirectionFromGlobalDirection();

// if local could not be updated from global
    if(!vf_direction || !inc_direction){
        v1->ReverseVFLocalDirection();
        v2->ReverseVFLocalDirection();
        v3->ReverseVFLocalDirection();
        if(v1->VertexOwnInclusion()){
            v1->GetInclusion()->Reverse_Direction();
        }
        if(v2->VertexOwnInclusion()){
            v2->GetInclusion()->Reverse_Direction();
        }
        if(v3->VertexOwnInclusion()){
            v3->GetInclusion()->Reverse_Direction();
        }
        
        if(ClosingAHole){
            KillATriangle(l2);
        }
        else{
            KillALink(newlink);
        }
        return false;
    }

    
    enew = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
    enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
    enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
    // vector field contributions
    enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v1);
    enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v2);
    enew += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v3);



    //=== inclusion interaction energy
    for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
    }
    for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
    }
    for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
    }

    if(!ClosingAHole){

        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);
        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, v1->m_pPrecedingEdgeLink);
        }
        // the new created link
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(newlink);
    
        for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
            enew +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, newlink);
        }
    }

    double diff_energy = enew - eold;
    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    if(exp(-m_Beta  *diff_energy + m_DBeta ) > thermal ) {
        
        m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
        return true;
    }
    else
    {
        if(ClosingAHole){
            KillATriangle(l2);
        }
        else{
            KillALink(newlink);
        }
        
        v1->ReverseVFLocalDirection();
        v2->ReverseVFLocalDirection();
        v3->ReverseVFLocalDirection();
        if(v1->VertexOwnInclusion()){
            v1->GetInclusion()->Reverse_Direction();
        }
        if(v2->VertexOwnInclusion()){
            v2->GetInclusion()->Reverse_Direction();
        }
        if(v3->VertexOwnInclusion()){
            v3->GetInclusion()->Reverse_Direction();
        }
        
        double e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
        e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
        e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
        e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v1);
        e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v2);
        e += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(v3);

        {
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it){
                e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it){
                e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer<m_No_VectorFields_Per_V; vf_layer++) {
                    (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it){
                e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                for( int vf_layer = 0; vf_layer< m_No_VectorFields_Per_V; vf_layer++) {
                    (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
        }
        e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);
        for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++) {
            (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, v1->m_pPrecedingEdgeLink);
        }
        return false;

    }
    
    return true;
}
// an atempt to create a link                       //        v3                      v3
// between v1 and v2                                // ml3  /   \  ml2---> ml3,l3  // T \\  l2, ml2   ml3(v1,v3)  ml2(v3,v2)
                                                  //     v1      v2               v1---->v2
links* OpenEdgeEvolutionWithConstantVertex::CreateALink(vertex *v1) {           //     l1
// l1 will be created and
    
    if(m_pGhostT.size()<1 && m_pGhostL.size()<3){
        std::cout<<" error 882---> we do not have enough storage for the new trinagle and links \n";
        exit(0);
    }
    
    links* ml3 = v1->GetEdgeLink();
    vertex *v3 = ml3->GetV2();
    links* ml2 = v3->GetEdgeLink();
    vertex *v2 = ml2->GetV2();
    links *l1 = m_pGhostL[0];
    links *l2 = m_pGhostL[1];
    links *l3 = m_pGhostL[2];
    triangle *tre = m_pGhostT[0];  // one trinagule is created

//---- create the triangle
    tre->UpdateVertex(v1,v2,v3);
    m_pGhostT.erase(m_pGhostT.begin());
    AddtoTriangleList(tre, m_pActiveT);
    l1->UpdateTriangle(tre);
    l2->UpdateTriangle(tre);
    l3->UpdateTriangle(tre);
    v1->AddtoTraingleList(tre);
    v2->AddtoTraingleList(tre);
    v3->AddtoTraingleList(tre);

//----- update the links
      //-- update l1 and l2 that have became a surface edge
    RemoveFromLinkList(ml2, m_pEdgeL);
    RemoveFromLinkList(ml3, m_pEdgeL);
    AddtoLinkList(ml2, m_pRightL);
    AddtoLinkList(ml3, m_pRightL);
    ml2->m_LinkType = 0;
    ml3->m_LinkType = 0;
    l2->m_LinkType = 0;
    l3->m_LinkType = 0;
    l1->m_LinkType = 1;
    ml2->UpdateMirrorFlag(true);
    ml3->UpdateMirrorFlag(true);
    l2->UpdateMirrorFlag(true);
    l3->UpdateMirrorFlag(true);
    l1->UpdateMirrorFlag(false);
    ml2->UpdateMirrorLink(l2);
    ml3->UpdateMirrorLink(l3);
    l2->UpdateMirrorLink(ml2);
    l3->UpdateMirrorLink(ml3);
    AddtoLinkList(l2, m_pLeftL);   // since we just added ml1 to m_pRightL
    AddtoLinkList(l3, m_pLeftL);
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());    // three times for l1,l2,l1
    AddtoLinkList(l1, m_pEdgeL);
    AddtoLinkList(l1, m_pActiveL);
    AddtoLinkList(l2, m_pActiveL);
    AddtoLinkList(l3, m_pActiveL);

    l1->UpdateNeighborLink1(l2);
    l1->UpdateNeighborLink2(l3);
    l2->UpdateNeighborLink1(l3);
    l2->UpdateNeighborLink2(l1);
    l3->UpdateNeighborLink1(l1);
    l3->UpdateNeighborLink2(l2);
    l1->UpdateV(v1,v2,v3);
    l2->UpdateV(v2,v3,v1);
    l3->UpdateV(v3,v1,v2);
//----update vertices
    //----update v3 that has become a surf
    RemoveFromVertexList(v3, m_pEdgeV);   // only this vertex becames a surf
    AddtoVertexList(v3, m_pSurfV);
    v3->m_VertexType = 0;                          // meaning it is not an edge vertex anymore
    //----update other vertices
    v1->m_pEdgeLink = l1;
    v2->m_pPrecedingEdgeLink = l1;
    v1->AddtoNeighbourVertex(v2);
    v2->AddtoNeighbourVertex(v1);
    v1->AddtoLinkList(l1);
    v2->AddtoLinkList(l2);
    v3->AddtoLinkList(l3);
    
 //----- now we should update the geomtry of the affected v, l, t
    tre->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
   // l2->UpdateNormal();   // normal of the links with mirror should be updated
   // l3->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    l1->UpdateEdgeVector(m_pBox);   // l1 is an edge link, we need only the edge vector and length
    l2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    l3->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(v3);  // v3 is now a surface vertex

    return l1;
}

//         v3                    v3
//  l3  // T1 \\ l2  -->  ml3  /   \  ml2
//      v1-l1- v2             v1    v2
bool OpenEdgeEvolutionWithConstantVertex::KillALink(links *l1)
{
        triangle *tri = l1->GetTriangle();
        vertex *v1 = l1->GetV1();
        vertex *v2 = l1->GetV2();
        vertex *v3 = l1->GetV3();
        links *l2 = l1->GetNeighborLink1();
        links *l3 = l2->GetNeighborLink1();
        links *ml2 = l2->GetMirrorLink();
        links *ml3 = l3->GetMirrorLink();
    
//--- remove the triangle
    RemoveFromTriangleList(tri, m_pActiveT);
    AddtoTriangleList(tri,m_pGhostT);
          //--- remove the trinagle from all the vertcies list
    v1->RemoveFromTraingleList(tri);
    v2->RemoveFromTraingleList(tri);
    v3->RemoveFromTraingleList(tri);
    
//---- remove the links l1, l2 and l3 ; note the mirors remain allive but will become an edge link
    RemoveFromLinkList(l1, m_pActiveL);  // too expensive
    RemoveFromLinkList(l2, m_pActiveL);  // too expensive
    RemoveFromLinkList(l3, m_pActiveL);  // too expensive
    RemoveFromLinkList(l2, m_pRightL);  // too expensive
    RemoveFromLinkList(l3, m_pRightL);  // too expensive
    RemoveFromLinkList(l2, m_pLeftL);  // too expensive
    RemoveFromLinkList(l3, m_pLeftL);  // too expensive
    RemoveFromLinkList(l1, m_pEdgeL);    // this link does not exist in m_pMHL and m_pHL
    AddtoLinkList(l1,m_pGhostL);
    AddtoLinkList(l2,m_pGhostL);
    AddtoLinkList(l3,m_pGhostL);
    v1->RemoveFromLinkList(l1);
    v2->RemoveFromLinkList(l2);
    v3->RemoveFromLinkList(l3);
    
//--- convert the links into edge links ml2 and ml3
    RemoveFromLinkList(ml2, m_pRightL);  // too expensive; I should find a better way
    RemoveFromLinkList(ml3, m_pRightL);  // too expensive
    RemoveFromLinkList(ml2, m_pLeftL);  // too expensive
    RemoveFromLinkList(ml3, m_pLeftL);  // too expensive
    ml2->UpdateMirrorFlag(false);
    ml3->UpdateMirrorFlag(false);
    ml2->m_LinkType = 1;
    ml3->m_LinkType = 1;
    AddtoLinkList(ml2, m_pEdgeL);    // adding the two mirror into the edge
    AddtoLinkList(ml3, m_pEdgeL);   // adding the two mirror into the edge

//--- convert the three vertices into edge vertex,
        //--- only v2 needs an update: v1 and v3 are already an edge vertex.
    RemoveFromVertexList(v3, m_pSurfV);  // too expensive
    AddtoVertexList(v3, m_pEdgeV);
    v3->m_VertexType = 1; //
       //-- v1 and v2 are no longer connected
    v1->RemoveFromNeighbourVertex(v2);
    v2->RemoveFromNeighbourVertex(v1);
       // now we updated the vertex edge links and preceding link; note: v3 edge links will not be changed and also v1  preceding link will not be changed
        v1->m_pEdgeLink = ml3;
        v3->m_pEdgeLink = ml2;
        v3->m_pPrecedingEdgeLink = ml3;
        v2->m_pPrecedingEdgeLink = ml2 ;
//--- now we should update the geomtry of the affected v, l
        ml2->UpdateEdgeVector(m_pBox);   // edge vector should be updated
        ml3->UpdateEdgeVector(m_pBox);   // edge vector should be updated
  
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v1);  // v1 is still an edge vertex
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v2);  // // v2 is still an edge vertex
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v3);  // v3 is now an edge vertex
    
    return true;
    
}
// an atempt to create a triangle                    //       v3              v3
// using v1                                          // l3   /   \ l2  -->  // T \\
                                                  //       v1---v2        v1 ====v2
                                                    //        l1
triangle* OpenEdgeEvolutionWithConstantVertex::CloseATriangleHole(vertex *v1)
{


    links* l1 = v1->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    links* l2 = v2->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l3 = v3->m_pEdgeLink;

        // we added the triangle
    if(m_pGhostT.size()<1 || m_pGhostL.size()<3){
        std::cout<<" error 882---> we do not have enough storage for the new trinagles or links \n";
        exit(0);
    }

//-- new triangle and links that need to be created.
    links *ml1 = m_pGhostL[0];
    links *ml2 = m_pGhostL[1];
    links *ml3 = m_pGhostL[2];
    triangle *outtriangle = m_pGhostT[0];
//-- build the new triangle
    m_pGhostT.erase(m_pGhostT.begin());
    AddtoTriangleList(outtriangle,m_pActiveT);
    // this trinagle is made of v1, v2 and v3 but aniwise
    outtriangle->UpdateVertex(v2,v1,v3);
    // add the trinagle to the vertcies list
    v1->AddtoTraingleList(outtriangle);
    v2->AddtoTraingleList(outtriangle);
    v3->AddtoTraingleList(outtriangle);
    // note, the exsisting links belong to different trinagles. ml1, ml2, ml3 making this trinagle
    ml1->UpdateTriangle(outtriangle);
    ml2->UpdateTriangle(outtriangle);
    ml3->UpdateTriangle(outtriangle);
//-- build the ml links
    //-- removing from ghost containers
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml1
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml2
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml3
    //--add ml links to active
    AddtoLinkList(ml1,m_pActiveL);
    AddtoLinkList(ml2,m_pActiveL);
    AddtoLinkList(ml3,m_pActiveL);
    //--- to mhl
    AddtoLinkList(ml1,m_pLeftL);  // note: this is true because l1 is added to hl
    AddtoLinkList(ml2,m_pLeftL);
    AddtoLinkList(ml3,m_pLeftL);
    //-- next link
    ml1->UpdateNeighborLink1(ml3);
    ml1->UpdateNeighborLink2(ml2);
    ml2->UpdateNeighborLink1(ml1);
    ml2->UpdateNeighborLink2(ml3);
    ml3->UpdateNeighborLink1(ml2);
    ml3->UpdateNeighborLink2(ml1);
    //-- update their vertex
    ml1->UpdateV(v2,v1,v3);
    ml2->UpdateV(v3,v2,v1);
    ml3->UpdateV(v1,v3,v2);
    //== adding the m links to vertices list
    v1->AddtoLinkList(ml3);
    v3->AddtoLinkList(ml2);
    v2->AddtoLinkList(ml1);
    ml1->m_LinkType = 0; //
    ml2->m_LinkType = 0; //
    ml3->m_LinkType = 0; //
    ml1->UpdateMirrorFlag(true);
    ml2->UpdateMirrorFlag(true);
    ml3->UpdateMirrorFlag(true);
    ml1->UpdateMirrorLink(l1);
    ml2->UpdateMirrorLink(l2);
    ml3->UpdateMirrorLink(l3);
//-- update l1,l2,l3
    l1->m_LinkType = 0; //
    l2->m_LinkType = 0; //
    l3->m_LinkType = 0; //
    l1->UpdateMirrorFlag(true);
    l2->UpdateMirrorFlag(true);
    l3->UpdateMirrorFlag(true);
    l1->UpdateMirrorLink(ml1);
    l2->UpdateMirrorLink(ml2);
    l3->UpdateMirrorLink(ml3);
    //  l1 and l2 to  m_pHL
    AddtoLinkList(l1,m_pRightL);  //note  ml1  ml2 and ml3 are in mhl
    AddtoLinkList(l2,m_pRightL);
    AddtoLinkList(l3,m_pRightL);
    RemoveFromLinkList(l1,m_pEdgeL);     // they are part of normal links now
    RemoveFromLinkList(l2,m_pEdgeL);
    RemoveFromLinkList(l3,m_pEdgeL);
//== uodate vertices
    v1->m_VertexType = 0; // meaning it is not an edge vertex anymore
    v2->m_VertexType = 0; // meaning it is not an edge vertex anymore
    v3->m_VertexType = 0; // meaning it is not an edge vertex anymore
    //-- removing from vedge containers
    RemoveFromVertexList(v1,m_pEdgeV);   // v1, v2 , v3 is no longer an edge vertex
    RemoveFromVertexList(v2,m_pEdgeV);
    RemoveFromVertexList(v3,m_pEdgeV);
    //--- adding to the containers
    AddtoVertexList(v1,m_pSurfV);
    AddtoVertexList(v2,m_pSurfV);
    AddtoVertexList(v3,m_pSurfV);

// updating the geomtry
      // triangle area and normal
    outtriangle->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
    ml1->UpdateNormal();   // normal of the links with mirror should be updated
    ml2->UpdateNormal();   // normal of the links with mirror should be updated
    ml3->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    ml1->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    ml2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    ml3->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(v3);  // v3 is now a surface vertex

    

    return outtriangle;
}
// the triangle should belong to l1 (l1 get killed at the end)
// ---------------------------                      //               v3               v3
// between v1 and v2                                //     l3     //  \\ l2 --->     /   \
                                                  //             v1 ====v2          v1---v2
                                                         //          l1

// this function kills target_triangle and the assosciated links l1,l2,l3. It send them to the ghost area
bool OpenEdgeEvolutionWithConstantVertex::KillATriangle(links *l1){
//=== note: the mirror will be killed not l1,
    // this are the objects that will be killed
    links *ml1 = l1->GetMirrorLink();
    links *ml3 = ml1->GetNeighborLink1();   //    in the reverse function ml1->UpdateNeighborLink1(ml3);
    links *ml2 = ml3->GetNeighborLink1();   //    in the reverse function ml3->UpdateNeighborLink1(ml2);
    triangle *target_triangle = ml1->GetTriangle();

    //-- other objects
    links *l2 = ml2->GetMirrorLink();
    links *l3 = ml3->GetMirrorLink();
    vertex *v1 = ml1->GetV2();
    vertex *v2 = ml1->GetV1();      //  in the reverse function ml2->UpdateV(v3,v2,v1);
    vertex *v3 = ml1->GetV3();

//-- remove the triangle, only from active t and verices. Note the links that have this triangle will die anyway
    AddtoTriangleList(target_triangle,m_pGhostT);
    RemoveFromTriangleList(target_triangle,m_pActiveT);
    v1->RemoveFromTraingleList(target_triangle);
    v2->RemoveFromTraingleList(target_triangle);
    v3->RemoveFromTraingleList(target_triangle);
    
//--- remove the links ml1, ml2, ml3 links
    RemoveFromLinkList(ml1,m_pActiveL);
    RemoveFromLinkList(ml2,m_pActiveL);
    RemoveFromLinkList(ml3,m_pActiveL);
    RemoveFromLinkList(ml1,m_pLeftL);
    RemoveFromLinkList(ml2,m_pLeftL);
    RemoveFromLinkList(ml3,m_pLeftL);
    RemoveFromLinkList(ml1,m_pRightL);  // we know that they were added to mhl but for more generality
    RemoveFromLinkList(ml2,m_pRightL);
    RemoveFromLinkList(ml3,m_pRightL);
    AddtoLinkList(ml1,m_pGhostL);
    AddtoLinkList(ml2,m_pGhostL);
    AddtoLinkList(ml3,m_pGhostL);
// update l1,l2 and l3 links
    RemoveFromLinkList(l1,m_pRightL);  // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l2,m_pRightL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l3,m_pRightL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l1,m_pLeftL); // we know they were added to phl but for generality
    RemoveFromLinkList(l2,m_pLeftL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l3,m_pLeftL); // this link also has became an edge links, should be removed form the list
    AddtoLinkList(l1,m_pEdgeL);    // so they became an edge link
    AddtoLinkList(l2,m_pEdgeL);
    AddtoLinkList(l3,m_pEdgeL);
    l1->UpdateMirrorFlag(false);
    l2->UpdateMirrorFlag(false);
    l3->UpdateMirrorFlag(false);
    l1->m_LinkType = 1;
    l2->m_LinkType = 1;
    l3->m_LinkType = 1;
//-- convert the vertices into edge
    RemoveFromVertexList(v1,m_pSurfV);
    RemoveFromVertexList(v2,m_pSurfV);
    RemoveFromVertexList(v3,m_pSurfV);
    AddtoVertexList(v1,m_pEdgeV);
    AddtoVertexList(v2,m_pEdgeV);
    AddtoVertexList(v3,m_pEdgeV);
    v1->m_VertexType = 1;
    v2->m_VertexType = 1;
    v3->m_VertexType = 1;
    
    v1->RemoveFromLinkList(ml3);    // in the counter function    v1->AddtoLinkList(ml3);
    v2->RemoveFromLinkList(ml1);    // in the counter function    v2->AddtoLinkList(ml1);
    v3->RemoveFromLinkList(ml2);    // in the counter function    v3->AddtoLinkList(ml2);

    v1->m_pEdgeLink = l1;    // for here is not needed (since v1 was an edge vector). But in a general case, v1 may not have an edge vector
    v2->m_pEdgeLink = l2;
    v3->m_pEdgeLink = l3;
    v1->m_pPrecedingEdgeLink = l3;
    v2->m_pPrecedingEdgeLink = l1;
    v3->m_pPrecedingEdgeLink = l2;
      // now we should update the geomtry of the affected v1,v2,v3, ml1,ml2,ml3,
    l1->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    l2->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    l3->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v3);  // v3 is now a surface vertex

    return true;
}
double  OpenEdgeEvolutionWithConstantVertex::SystemEnergy() {
    double en = 0;
    /*
    std::vector<vertex *> ActiveV =  m_pSurfV;
    std::vector<triangle *> pActiveT =  m_pActiveT;
    std::vector<links *> mLink =  m_pHL;
    std::vector<links *>  pEdgeL =  m_pEdgeL;
    std::vector<vertex *> EdgeV  =  m_pEdgeV;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = ( m_pHL).begin() ; it != ( m_pHL).end(); ++it){
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    for (std::vector<vertex *>::iterator it = ( m_pSurfV).begin() ; it != ( m_pSurfV).end(); ++it)
        (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(*it);
    
    en = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    */
    return en;
}
bool OpenEdgeEvolutionWithConstantVertex::Linkisvalid(vertex *v1) {
    // Check if the new link length is within the allowed range and if the angle of the new triangle
    // is acceptable with respect to the two other triangles.
    //           va -- v3----vb
    //            \n2 /n1\n3/
    //             v1--- v2
    //                DP
    
    // check if the new link length is within the allowed range and also if the angle of the new trinagule is fine with respect to the two other trinagules
    links* l2   = v1->GetEdgeLink();
    vertex *v3  = l2->GetV2();
    links* l1   = v3->GetEdgeLink();
    vertex *v2  = l1->GetV2();


    Vec3D DP = v2->GetPos() - v1->GetPos();
    
    for (int i=0;i<3;i++)
    if(fabs(DP(i)) > (*m_pBox)(i)/2.0)
    {
        if(DP(i) < 0)
            DP(i) = (*m_pBox)(i)+DP(i);
        else if(DP(i) > 0)
            DP(i) = DP(i)-(*m_pBox)(i);
    }
    double dist2 = DP.dot(DP,DP);
   // std::cout<<dist2<<"\n";
    if(dist2 < m_MinLength2 || dist2 > m_MaxLength2){
        return false;
    }
    
    triangle t(0,v1,v2,v3);
    t.UpdateNormal_Area(m_pBox);
    Vec3D n1 = t.GetNormalVector();
    Vec3D n2 = l2->GetTriangle()->GetNormalVector();
    Vec3D n3 = l1->GetTriangle()->GetNormalVector();
    
    if( n1.dot(n1,n2) < m_MinAngle){
        return false;
    }
    if( n1.dot(n1,n3) < m_MinAngle){
        return false;
    }
    
    return true;
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
// this can be written as template
void OpenEdgeEvolutionWithConstantVertex::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
std::string OpenEdgeEvolutionWithConstantVertex::CurrentState(){
    
    std::string state = AbstractOpenEdgeEvolution::GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + +" "+ Nfunction::D2S(m_Period) +" "+ Nfunction::D2S(m_NumberOfMovePerStep);

    return state;
}

