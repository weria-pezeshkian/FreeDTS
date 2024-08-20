#include <chrono>
#include <unordered_set>
#include <utility>
#include <ctime>
#include "Three_Edge_Scission.h"
#include "State.h"
#include "MESH.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */

Three_Edge_Scission::Three_Edge_Scission(int period, State *pState) :
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
                m_Beta(pState->GetSimulation()->GetBeta()),
                m_DBeta(pState->GetSimulation()->GetDBeta()),
                m_MinLength2(pState->GetSimulation()->GetMinL2()),
                m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
                m_MinAngle(pState->GetSimulation()->GetMinAngle()),
                m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex()),
                m_Box(pState->GetMesh()->Link2ReferenceBox())
{

}
Three_Edge_Scission::~Three_Edge_Scission(){
    
}
void Three_Edge_Scission::Initialize() {

    std::cout<<"---> note: Three_Edge Neck will be used to change the surface topology "<<std::endl;
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;
    m_Surface_Genus = 1-(m_pSurfV.size()-m_pLeftL.size()+m_pActiveT.size())/2;
}
bool Three_Edge_Scission::MCMove(int step) {
    
    if(m_Period == 0 ||  step%m_Period != 0){
        return false;
    }
    
     std::vector<pair_pot_triangle> pair_list  = FindNecks();
   // std::cout<<pair_list.size()<<" number of pot trinagles \n";
   if(pair_list.size() != 0) {// ScissionByMC
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(pair_list.size());
        pair_pot_triangle pair_T = pair_list[n];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(ScissionByMC(pair_T, thermal)){
            m_AcceptedMoves++;
        }
    } //  if(pair_list.size() != 0) end ScissionByMC
    std::clock_t start = std::clock();
    std::vector<fussion_site> all_possible_sites = FindPotentialFussionSites();
    std::clock_t end = std::clock();
    // Calculate the duration
       double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
       // Print the time taken
       std::cout << "Time taken to execute FindPotentialFussionSites: " << duration << " seconds " << all_possible_sites.size() << std::endl;
    if(all_possible_sites.size() != 0) {// FussionByMove
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(pair_list.size());
        fussion_site pair_T = all_possible_sites[n];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(FussionByMove(pair_T, thermal)){
            m_AcceptedMoves++;
        }
    }//     if(all_possible_sites.size() != 0) {// end FussionByMove

    m_Surface_Genus = 1 - (m_pSurfV.size()-m_pLeftL.size()+m_pActiveT.size())/2;
    return true;
}
bool Three_Edge_Scission::FussionByMove(fussion_site &pair_tri, double thermal){
    return false;
}
std::vector<fussion_site> Three_Edge_Scission::FindPotentialFussionSites(){
    
    std::vector<fussion_site> all_possible_sites;
    
    //MESH::SquareDistanceBetweenTwoVertices(p_v1, p_v2, m_Box);

   // for (std::vector<links *>::iterator it = m_pRightL.begin() ; it != m_pRightL.end(); ++it){

        for(int i = 0; i<m_pActiveL.size(); i++){
            for(int j = 0; j<m_pActiveL.size(); j++){
               // if(all_possible_sites.size())
                 //   break;
                
                links *l1 = m_pActiveL[i];
                links *l2 = m_pActiveL[j];
                //---- check if they are close
                fussion_site p_T;
                if(CheapScane(l1,l2, p_T)){
                    continue;
                }

                if(BuildScane(p_T)){
                    all_possible_sites.push_back(p_T);
                    
                    l1->GetV1()->UpdateGroup(1);
                    l1->GetV2()->UpdateGroup(1);
                    l1->GetV3()->UpdateGroup(1);
                    l2->GetV1()->UpdateGroup(2);
                    l2->GetV2()->UpdateGroup(2);
                    l2->GetV3()->UpdateGroup(2);

                }
            //====
            }
        }

    
    
    //n1.n2<0??
    
    return all_possible_sites;
}
bool Three_Edge_Scission::BuildScane(fussion_site &p_T) {

    
    
    return true;
}
bool Three_Edge_Scission::CheapScane(links *l1, links *l2, fussion_site &p_T) {
    /**
     * @brief Determines if two links are close to each other based on their normal vectors and vertex proximity.
     *
     * This function checks if two given links are close by first comparing the normal vectors of the triangles
     * they belong to (because fussion trinagles should have this characters). If the dot product of the normal vectors is positive, the links are considered close.
     * Additionally, it checks the proximity of the vertices of the two links. If any vertex of the second link
     * is a neighbor of any vertex of the first link, the links are considered close.
     *
     * @param l1 Pointer to the first link.
     * @param l2 Pointer to the second link.
     * @return True if the links are considered close, false otherwise.
     */
#if DEVELOPMENT_MODE == Enabled
    std::cout << "DEVELOPMENT_MODE ID 00345: This function can be made better\n";
#endif

    // If the dot product of the normal vectors is positive, return true
    if (Vec3D::dot(l1->GetTriangle()->GetNormalVector(), l2->GetTriangle()->GetNormalVector()) > 0) {
        return true;
    }

    double dist_33 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV3(), l2->GetV3(), m_Box);
    if(dist_33 > m_MaxLength2){
        return true;
    }
    
    if(l1->GetV1() == l2->GetV1() ||
       l1->GetV1() == l2->GetV2() ||
       l1->GetV1() == l2->GetV3() ){
        return true;
    }
    if(l1->GetV2() == l2->GetV1() ||
       l1->GetV2() == l2->GetV2() ||
       l1->GetV2() == l2->GetV3() ){
        return true;
    }
    if(l1->GetV3() == l2->GetV1() ||
       l1->GetV3() == l2->GetV2() ||
       l1->GetV3() == l2->GetV3() ){
        return true;
    }
    if(l1->GetV1() == l2->GetMirrorLink()->GetV3() ||
       l1->GetV2() == l2->GetMirrorLink()->GetV3() ||
       l1->GetV3() == l2->GetMirrorLink()->GetV3() ){
        return true;
    }
    if(l1->GetV1() == l2->GetNeighborLink1()->GetMirrorLink()->GetV3() ||
       l1->GetV2() == l2->GetNeighborLink1()->GetMirrorLink()->GetV3() ||
       l1->GetV3() == l2->GetNeighborLink1()->GetMirrorLink()->GetV3() ){
        return true;
    }
    if(l1->GetV1() == l2->GetNeighborLink2()->GetMirrorLink()->GetV3() ||
       l1->GetV2() == l2->GetNeighborLink2()->GetMirrorLink()->GetV3() ||
       l1->GetV3() == l2->GetNeighborLink2()->GetMirrorLink()->GetV3() ){
        return true;
    }
    if(l2->GetV1() == l1->GetMirrorLink()->GetV3() ||
       l2->GetV2() == l1->GetMirrorLink()->GetV3() ||
       l2->GetV3() == l1->GetMirrorLink()->GetV3() ){
        return true;
    }
    if(l2->GetV1() == l1->GetNeighborLink1()->GetMirrorLink()->GetV3() ||
       l2->GetV2() == l1->GetNeighborLink1()->GetMirrorLink()->GetV3() ||
       l2->GetV3() == l1->GetNeighborLink1()->GetMirrorLink()->GetV3() ){
        return true;
    }
    if(l2->GetV1() == l1->GetNeighborLink2()->GetMirrorLink()->GetV3() ||
       l2->GetV2() == l1->GetNeighborLink2()->GetMirrorLink()->GetV3() ||
       l2->GetV3() == l1->GetNeighborLink2()->GetMirrorLink()->GetV3() ){
        return true;
    }


    double dist_12 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV1(), l2->GetV2(), m_Box);
    double dist_21 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV2(), l2->GetV1(), m_Box);

    if(dist_21 > m_MaxLength2 || dist_12 > m_MaxLength2){
        return true;
    }
    
    double dist_11 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV1(), l2->GetV1(), m_Box);
    double dist_22 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV2(), l2->GetV2(), m_Box);
    if( dist_11> m_MaxLength2 && dist_22> m_MaxLength2){
        return true;
    }
    
    double dist_13 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV1(), l2->GetV3(), m_Box);
    double dist_32 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV3(), l2->GetV2(), m_Box);
    if (dist_13 > m_MaxLength2 && dist_32> m_MaxLength2 ){
            return true;
    }

    double dist_23 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV2(), l2->GetV3(), m_Box);
    double dist_31 = MESH::SquareDistanceBetweenTwoVertices(l1->GetV3(), l2->GetV1(), m_Box);

    if (dist_23> m_MaxLength2 && dist_31> m_MaxLength2){
            return true;
    }
    p_T.l1 = l1;
    p_T.l2 = l2;
    p_T.dist[0][0] = dist_11;
    p_T.dist[0][1] = dist_12;
    p_T.dist[0][2] = dist_13;
    p_T.dist[1][0] = dist_21;
    p_T.dist[1][1] = dist_22;
    p_T.dist[1][2] = dist_23;
    p_T.dist[2][1] = dist_32;
    p_T.dist[2][2] = dist_33;
    p_T.dist[2][2] = dist_33;

    return false;
}
bool Three_Edge_Scission::ScissionByMC(pair_pot_triangle &pair_t, double thermal){
    /**
     * @brief Perform a scission operation on a neck by getting a potential pair of triangles and determine its acceptance based on Metropolis criteria.
     *
     * @param pair_t A potential pair of triangles that will be created after  the scission.
     * @param thermal The thermal factor used in the Metropolis acceptance criterion.
     * @return true if the scission is accepted, false otherwise.
     *
     * This function calculates the energy before and after performing a scission operation on a pair of triangles.
     * It then uses the Metropolis criterion to decide whether to accept the new configuration or revert to the old one.
     */
    // Check if there are enough links and triangles in the repository
    if (m_pGhostT.size() < 4 || m_pGhostL.size() < 4) {
        std::cout << " --->note: the number of the links and triangles in the repository is not enough, restart the simulations \n";
        return false;
    }
    
    double new_energy = 0;
    double old_energy = 0;
    
    //---> get all the links
        vertex *v11 = pair_t.PT1.pv1;
        vertex *v12 = pair_t.PT1.pv2;
        vertex *v13 = pair_t.PT1.pv3;
        vertex *v21 = pair_t.PT2.pv1;
        vertex *v22 = pair_t.PT2.pv2;
        vertex *v23 = pair_t.PT2.pv3;
    
        v11->EnergyCopy();
        v12->EnergyCopy();
        v13->EnergyCopy();
        v21->EnergyCopy();
        v22->EnergyCopy();
        v23->EnergyCopy();
//---> calculate old energies
    //---- effected vertex energy
    old_energy += v11->GetEnergy();
    old_energy += v12->GetEnergy();
    old_energy += v13->GetEnergy();
    old_energy += v21->GetEnergy();
    old_energy += v22->GetEnergy();
    old_energy += v23->GetEnergy();
    
//--> get the links that may be effected by the cut
    std::vector <links *> Clinks = pair_t.ConnectingLinks;
    std::vector <triangle *> Ctriangles = pair_t.ConnectingTriangles;
    
    // find the links in which there interaction energy changes
    std::vector<links*> Affected_links_old = GetEdgesWithInteractionChange(pair_t);
    for (std::vector<links *>::iterator it = Affected_links_old.begin() ; it != Affected_links_old.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }
//----> Perform the scission
    std::vector<triangle *> pair_tri = DoAScission(pair_t);

//----> Calculate new energy
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v11);
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v12);
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v13);
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v21);
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v22);
    new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v23);

        std::vector<links*> Affected_links_new = GetEdgesWithInteractionChange(pair_t);
        for (std::vector<links *>::iterator it = Affected_links_new.begin() ; it != Affected_links_new.end(); ++it){
            new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        
            if(pair_t.PT1.pv1->GetNumberOfVF() != 0 ){
                for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
        }
    
    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy;
    double U = m_Beta * tot_diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > thermal ) {
        //--- Accepted
        m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
        /*
         
            // area and voulme effects

        }*/
        return true;
    }
    else{ // reject the move
        
        // Rejected, reverse the scission
        ReverseAScission(pair_t, pair_tri[0], pair_tri[1]);
        
        // Recalculate energies for consistency (although not needed for return value)
        v11->ReverseEnergyCopy();
        v12->ReverseEnergyCopy();
        v13->ReverseEnergyCopy();
        v21->ReverseEnergyCopy();
        v22->ReverseEnergyCopy();
        v23->ReverseEnergyCopy();

        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        //std::vector<links*> Affected_links_old = GetEdgesWithInteractionChange(pair_t);
        for (std::vector<links *>::iterator it = Affected_links_old.begin() ; it != Affected_links_old.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
            (*it)->Reverse_VFInteractionEnergy();
        }


        return false;
    }
    
    return true;
}
// it creates a triangle and place it to the active trinagles list
triangle * Three_Edge_Scission::CreateATriangleFromAPotentialTriangle(pot_triangle &p1) {
//---> initate the copy of the three links for the potential triangle
    (p1.pl1)->SetCopy();
    (p1.pl2)->SetCopy();
    (p1.pl3)->SetCopy();
    
//---> create the triangle
    //--- select the last ghost triangle
    triangle *gt1 = m_pGhostT[m_pGhostT.size()-1];
    //--- update its vertex
    gt1->UpdateVertex(p1.pv1,p1.pv2,p1.pv3);

    //--- put the triangle from ghost to active
    m_pActiveT.push_back(gt1);
    m_pGhostT.pop_back();
    
//---> now make the changes in the edges
    //--- update their triangles
    (p1.pl1)->UpdateTriangle(gt1);
    (p1.pl2)->UpdateTriangle(gt1);
    (p1.pl3)->UpdateTriangle(gt1);
    //---- update their next two edges
    p1.pl1->UpdateNeighborLink1(p1.pl2);
    p1.pl1->UpdateNeighborLink2(p1.pl3);
    p1.pl2->UpdateNeighborLink1(p1.pl3);
    p1.pl2->UpdateNeighborLink2(p1.pl1);
    p1.pl3->UpdateNeighborLink1(p1.pl1);
    p1.pl3->UpdateNeighborLink2(p1.pl2);
    //--- update their v3 (v2 and v1 remain the same)
    p1.pl1->UpdateV3(p1.pv3);
    p1.pl2->UpdateV3(p1.pv1);
    p1.pl3->UpdateV3(p1.pv2);
    
//---> update the vertices triangle list
    p1.pv1->AddtoTraingleList(gt1);
    p1.pv2->AddtoTraingleList(gt1);
    p1.pv3->AddtoTraingleList(gt1);
    
//---> update geometry; we cannot update the geometry of the the vertices because they have some edges that should be removed
    //--- update the normal and area of the triangle
    gt1->UpdateNormal_Area(&m_Box);   //trinagule normal and area should be obtained
    //--- update the three links notmal and shape Operator
    (p1.pl1)->UpdateNormal();
    (p1.pl2)->UpdateNormal();
    (p1.pl3)->UpdateNormal();
    (p1.pl1)->UpdateShapeOperator(&m_Box);
    (p1.pl2)->UpdateShapeOperator(&m_Box);
    (p1.pl3)->UpdateShapeOperator(&m_Box);

    return gt1;
}
// this function cuts the neck made of p1 and p2, i.e., pair
std::vector <triangle *> Three_Edge_Scission::DoAScission(pair_pot_triangle &pair){
    
    std::vector <triangle *> createdtriangles;
    if(m_pGhostT.size()<2){
        std::cout<<" ---> [not an error] not enough reserved trinagles; you may restart the simulation "<<std::endl;
        exit(0);
    }
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
//---> generate the two triangles and store them in the return vector
    //--- create the triangles
    createdtriangles.push_back(CreateATriangleFromAPotentialTriangle(p1));
    createdtriangles.push_back(CreateATriangleFromAPotentialTriangle(p2));

//---> remove the links that connect the two potential trinagles; we just remove them from different lists and send them to ghost
    //-- these edges are stored in the Clinks that are obtained in when the pair is created
    //-- note, the mirror links do not exist in this list, also their triangles should be deleted
    //-- these two triangles now have been created
    std::vector <links *> Clinks = pair.ConnectingLinks;       // this does not include the mirror links
    for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){
            MakeALinkGhost(*it);
    }
    std::vector <triangle *> all_triangle = pair.ConnectingTriangles;
    for (std::vector<triangle*>::iterator it = all_triangle.begin() ; it != all_triangle.end(); it++){
        MakeATriangleGhost(*it);
    }
//--> update geometry
    //-- update geometry of the 6 vertices
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv3);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv3);

    return createdtriangles;
}
// this is the exact reverse action of DoAScission; different from DoAFussion
bool Three_Edge_Scission::ReverseAScission(pair_pot_triangle &pair , triangle *t1, triangle *t2){

//--- getting p1 and p2
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
//-- we take the triangle3 t1 and t2 to the ghost and remove them from associated vertices and edges
    //-- send t1 and t2 to the ghost
    RemoveFromTriangleList(t1, m_pActiveT);
    RemoveFromTriangleList(t2, m_pActiveT);
    m_pGhostT.push_back(t1);
    m_pGhostT.push_back(t2);
    //--- remove t1 and t2 from their v->t list
    (t1->GetV1())->RemoveFromTraingleList(t1);
    (t1->GetV2())->RemoveFromTraingleList(t1);
    (t1->GetV3())->RemoveFromTraingleList(t1);
    (t2->GetV1())->RemoveFromTraingleList(t2);
    (t2->GetV2())->RemoveFromTraingleList(t2);
    (t2->GetV3())->RemoveFromTraingleList(t2);
    //--- reverse the edges to previous value
    (p1.pl1)->Reverse2PreviousCopy();
    (p1.pl2)->Reverse2PreviousCopy();
    (p1.pl3)->Reverse2PreviousCopy();
    (p2.pl1)->Reverse2PreviousCopy();
    (p2.pl2)->Reverse2PreviousCopy();
    (p2.pl3)->Reverse2PreviousCopy();

     //
//==== from the pair, we have all the links that were used to connect t1 and t2, we can recover them well
    std::vector <links *> Clinks = pair.ConnectingLinks;       // this does not include the mirror links
    for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){
        ActivateAGhostLink(*it);
    }
    std::vector <triangle *> all_triangle = pair.ConnectingTriangles;
    for (std::vector<triangle*>::iterator it = all_triangle.begin() ; it != all_triangle.end(); it++){
        ActivateAGhostTriangle(*it);
    }
//--> update geometry
    //-- update geometry of the 6 vertices
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv3);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv3);
    
    return true;
}
void Three_Edge_Scission::MakeALinkGhost(links *p_links){
    /**
     * @brief Make a link a "ghost" link, effectively deactivating it while preserving its information for potential reversal.
     *
     * This function deactivates a given link and its mirror by moving them from the active lists to the ghost list.
     * The link and its associated triangles are removed from various lists including the link list of vertices and
     * triangle lists. This process can be reversed without harm, making the link and its information valid if needed.
     *
     * @param p_links A pointer to the link to be made a ghost.
     *
     * @note The function performs checks in debug mode to ensure the link is not already in the ghost list. If the link
     * or its mirror link is found in the ghost list, an error message is printed and the function exits.
     *
     * @warning Ensure that the link and its mirror link are properly initialized and not already deactivated before
     * calling this function.
     */
    
#if DEBUG_MODE == Enabled
    // Check if the link or its mirror is already a ghost
    if (std::find(m_pGhostL.begin(), m_pGhostL.end(), p_links) != m_pGhostL.end() ||
        std::find(m_pGhostL.begin(), m_pGhostL.end(), p_links->GetMirrorLink()) != m_pGhostL.end()) {
        std::cout << "---> error 83732: such a function should have been called for this edge \n";
        return;
    }
#endif
    
    links *p_mlinks = p_links->GetMirrorLink();
    //-- remove the edge and its mirror
    RemoveFromLinkList(p_links,  m_pActiveL);
    RemoveFromLinkList(p_links,  m_pRightL);
    RemoveFromLinkList(p_links,  m_pLeftL);
    RemoveFromLinkList(p_mlinks,  m_pActiveL);
    RemoveFromLinkList(p_mlinks,  m_pRightL);
    RemoveFromLinkList(p_mlinks,  m_pLeftL);
    
    // Add the link and its mirror to the ghost list
    m_pGhostL.push_back(p_links);
    m_pGhostL.push_back(p_mlinks);
    
    //-- remove the edge from the list of the vertices
    (p_links->GetV1())->RemoveFromLinkList(p_links);
    (p_links->GetV2())->RemoveFromLinkList(p_mlinks);

    
    //-- make sure that  v1 and v2 are not connected by v-list; we do not need this for mirror link
    (p_links->GetV1())->RemoveFromNeighbourVertex(p_links->GetV2());
    (p_links->GetV2())->RemoveFromNeighbourVertex(p_links->GetV1());

    // what about ghost link
    return;
}
void Three_Edge_Scission::MakeATriangleGhost(triangle *p_tri){

    //-- remove the assosciated triangle
    RemoveFromTriangleList(p_tri, m_pActiveT);
    m_pGhostT.push_back(p_tri);
    //-- remove the associated triangle from the vertices list
    p_tri->GetV1()->RemoveFromTraingleList(p_tri);
    p_tri->GetV2()->RemoveFromTraingleList(p_tri);
    p_tri->GetV3()->RemoveFromTraingleList(p_tri);
    
    // what about ghost link
    return;
}
void Three_Edge_Scission::ActivateAGhostTriangle(triangle *p_tri){

    RemoveFromTriangleList(p_tri,m_pGhostT);
    AddtoVectorCarefully(p_tri, m_pActiveT);
    //--- add the triangles to the vertex list
    p_tri->GetV1()->AddtoTriangleListCarefully(p_tri);
    p_tri->GetV2()->AddtoTriangleListCarefully(p_tri);
    p_tri->GetV3()->AddtoTriangleListCarefully(p_tri);

    return;
}
void Three_Edge_Scission::ActivateAGhostLink(links *p_links){
    /**
     * @brief Reactivate a ghost link, restoring its activity and associated structures.
     *
     * This function reactivates a previously ghosted link and its mirror by moving them from the ghost list to the active lists.
     * It restores the link and its associated triangles to various active lists, ensuring they are properly reconnected within
     * the mesh or network structure. The function updates the link lists of the vertices and reestablishes the neighbor relationships between the vertices connected by the link.
     *
     * @param p_links A pointer to the link to be reactivated.
     *
     * @note The function assumes the link and its mirror link are currently in the ghost list. It carefully adds the
     * link and its associated triangle back to active lists, and ensures all necessary connections are reestablished.
     *
     * @warning Ensure that the link and its mirror link are correctly initialized and currently ghosted before calling
     * this function. Misuse may lead to inconsistencies in the data structure.
     */
    
    
    links *p_mlinks = p_links->GetMirrorLink();

    //--removing it and its mirror from ghost and add to real containors
    m_pActiveL.push_back(p_links);
    m_pActiveL.push_back(p_mlinks);
    m_pLeftL.push_back(p_links);
    m_pRightL.push_back(p_mlinks);

    RemoveFromLinkList(p_links,m_pGhostL);
    RemoveFromLinkList(p_mlinks,m_pGhostL);

    //-- add the links to the vertex linklist
    (p_links->GetV1())->AddtoLinkListCarefully(p_links);
    (p_links->GetV2())->AddtoLinkListCarefully(p_links->GetMirrorLink());
    //--- add the triangles to the vertex list
    //--- make v1 and v2 nighbour again
    (p_links->GetV1())->AddtoNeighbourVertexCarefully(p_links->GetV2());
    (p_links->GetV2())->AddtoNeighbourVertexCarefully(p_links->GetV1());
    
    
    return;
}
//== this function get a mesh and search through it and finds possible fission sites
std::vector<pair_pot_triangle> Three_Edge_Scission::FindNecks(){

//---- this part searches through all the links and finds possible triangles that do not yet exsist. This means, it finds triple vertices (v1,v2,v3) that are connected by edges but such trinagle is not defined.
    
    int id = 0;
    std::vector<pot_triangle> list;
    std::vector<links*>  all_link =  m_pActiveL;
    for (std::vector<links*>::iterator il1 = all_link.begin() ; il1 != all_link.end(); il1++){
        vertex *pv1 = (*il1)->GetV1();
        vertex *pv2 = (*il1)->GetV2();
        std::vector<links*>  v2_l = pv2->GetVLinkList();
        for (std::vector<links*>::iterator il2 = v2_l.begin() ; il2 != v2_l.end(); il2++){
            vertex *pv3 = (*il2)->GetV2();
            std::vector<links*>  v3_l = pv3->GetVLinkList();
                for (std::vector<links*>::iterator il3 = v3_l.begin() ; il3 != v3_l.end(); il3++){
                    if(pv1==(*il3)->GetV2() && (*il3)->GetV3()!=pv2 && ((*il3)->GetMirrorLink())->GetV3()!=pv2 ){
                        if( pv2->m_VertexType==0 && pv3->m_VertexType==0 && pv1->m_VertexType==0){
                            // pv1, pv2, pv3 are our vertices
                                                    
                            pot_triangle PotT;
                            PotT.id = id; id++;
                            PotT.cid = -1;
                            PotT.pv1 = pv1; PotT.pv2 = pv2; PotT.pv3 = pv3;
                            PotT.pl1= (*il1); PotT.pl2= (*il2); PotT.pl3= (*il3);
                            list.push_back(PotT);
                            //--- we remove these links from all link containor so we do not search through them again
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il1)->GetMirrorLink()), all_link.end());

                        }
                    }
            }
        }
    }
//===================
std::vector <pair_pot_triangle> pair_list;
    int idpair=0;
for (int i=0;i<list.size();i++)
{
    for (int j=i+1;j<list.size();j++){
        pair_pot_triangle temp;
        if(Is_A_Neck(list[i],list[j], temp)){
            temp.id = idpair;
            idpair++;
            pair_list.push_back(temp);
        }
    }
}
    
    return pair_list;
}
//-- connected_2pot_triangles function check if the two potential trinagles are well connected for a fission. not, T1 and T2 do not exist yet, but they can apear if v1,v2,v3 get
//          v1------------v4            disconnected from v4, v5, v6. This function checks for such cases
//         /T1\         / T2\
//        v2--v3-------v4---v6
bool Three_Edge_Scission::Is_A_Neck(pot_triangle potT1, pot_triangle potT2, pair_pot_triangle & neck) {
    /**
     * @brief Checks if two potential triangles (pot_triangle) are well-connected to form a neck for fission.
     *
     * This function determines if two potential triangles (pot_triangle) form a neck for fission. The triangles do not yet exist
     * but can appear if vertices v1, v2, and v3 are disconnected from v4, v5, and v6. The function ensures that the vertices of the first
     * potential triangle are connected to at least one vertex of the second potential triangle and vice versa.
     * If the triangles are well-connected, it updates their connection IDs and returns a pair of potential triangles
     * with connecting links and triangles.
     *
     * The function performs the following steps:
     * 1. Checks if either triangle is already connected to another potential triangle.
     * 2. Ensures the two potential triangles do not share any vertices.
     * 3. Verifies if the vertices of the first triangle are connected to the vertices of the second triangle.
     * 4. Checks the orientation of the triangles to ensure valid neck formation.
     * 5. Finds all links that connect the vertices of the first triangle to the second triangle and updates the connection information.
     *
     * @param potT1 First potential triangle.
     * @param potT2 Second potential triangle.
     * @param neck Reference to the pair_pot_triangle struct to store the result.
     * @return true if the potential triangles form a neck, false otherwise.
     //          v1------------v4            disconnected from v4, v5, v6. This function checks for such cases
     //         /T1\              / T2\
     //        v2--v3-------v4---v6
     */
    neck.state = false;
    
//--- Check if either triangle is already connected to another potential triangle.
    // Maybe this is not needed as we only make one fission
    if(potT2.cid != -1 || potT1.cid != -1){
        return false;
    }
    
//---> if the two potential trinagles share a vertex, then they should be discarded
    // Ensure the two potential triangles do not share any vertices.
    if (potT2.pv1 == potT1.pv1 || potT2.pv1 == potT1.pv2 || potT2.pv1 == potT1.pv3 ||
        potT2.pv2 == potT1.pv1 || potT2.pv2 == potT1.pv2 || potT2.pv2 == potT1.pv3 ||
        potT2.pv3 == potT1.pv1 || potT2.pv3 == potT1.pv2 || potT2.pv3 == potT1.pv3) {
        return false;
    }

//----> this check if the potential triangles are connected or not
    links* linksToCheck1[] = { potT1.pl1, potT1.pl2, potT1.pl3 };
    for (int i = 0; i < 3; ++i) {
        links* oe_l = linksToCheck1[i];
        links* me_l = oe_l->GetMirrorLink();

        vertex* v3_oe = oe_l->GetV3();
        vertex* v3_me = me_l->GetV3();

        if (v3_oe != potT2.pv1 && v3_oe != potT2.pv2 && v3_oe != potT2.pv3 &&
            v3_me != potT2.pv1 && v3_me != potT2.pv2 && v3_me != potT2.pv3) {
            return false;
        }
    }
    links* linksToCheck2[] = { potT2.pl1, potT2.pl2, potT2.pl3 };
    for (int i = 0; i < 3; ++i) {
        links* oe_l = linksToCheck2[i];
        links* me_l = oe_l->GetMirrorLink();

        vertex* v3_oe = oe_l->GetV3();
        vertex* v3_me = me_l->GetV3();

        if (v3_oe != potT1.pv1 && v3_oe != potT1.pv2 && v3_oe != potT1.pv3 &&
            v3_me != potT1.pv1 && v3_me != potT1.pv2 && v3_me != potT1.pv3) {
            return false;
        }
    }

    //--- here we should find all the links that connects the  vertices of P1 to P2.
    // how to return this
    if(CorrectOrientation(potT1,potT2) && CorrectOrientation(potT2,potT1)){ // this means that the p1 edges should be so that the removal edge disconect
                                                                 //vertices of p2 from p1, not any other vertices
        
        
         
        neck.PT1 = potT1;
        neck.PT2 = potT2;
        
        if(!CheckFaceAngle(neck)){
            //--> before doing more, lets check if this fission happens, could the angle be OK?
            return false;
        }
        std::vector <links *> Clinks;
        std::unordered_set <triangle *> Ctriangles;
        neck.state = true;

        // Iterate through the links connected to the first vertex of potT1
        std::vector<links*> n_p1links = (potT1.pv1)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p1links.begin(); it != n_p1links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.insert((*it)->GetTriangle());
            }
        }

        // Iterate through the links connected to the second vertex of potT1
        std::vector<links*> n_p2links = (potT1.pv2)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p2links.begin(); it != n_p2links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.insert((*it)->GetTriangle());
            }
        }

        // Iterate through the links connected to the third vertex of potT1
        std::vector<links*> n_p3links = (potT1.pv3)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p3links.begin(); it != n_p3links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.insert((*it)->GetTriangle());
            }
        }
        
        neck.ConnectingLinks = Clinks;
        for (std::unordered_set<triangle *>::iterator it = Ctriangles.begin(); it != Ctriangles.end(); ++it) {
            (neck.ConnectingTriangles).push_back(*it);
        }
    }
    else {
        std::cout<<"---> error2922: this should not happen \n";
    }
    
    return true;
}
bool Three_Edge_Scission::CheckFaceAngle(pair_pot_triangle &pair){
    
    triangle t1(0, pair.PT1.pv1, pair.PT1.pv2, pair.PT1.pv3);
    Vec3D n1 = t1.CalculateNormal(m_Box);
    Vec3D n1_1 = (pair.PT1.pl1)->GetMirrorLink()->GetTriangle()->GetNormalVector();
    Vec3D n1_2 = (pair.PT1.pl2)->GetMirrorLink()->GetTriangle()->GetNormalVector();
    Vec3D n1_3 = (pair.PT1.pl3)->GetMirrorLink()->GetTriangle()->GetNormalVector();

    if( n1.dot(n1,n1_1) < m_MinAngle ||
        n1.dot(n1,n1_2) < m_MinAngle ||
        n1.dot(n1,n1_3) < m_MinAngle ){
        return false;
    }
    
    triangle t2(0, pair.PT2.pv1, pair.PT2.pv2, pair.PT2.pv3);
    Vec3D n2 = t2.CalculateNormal(m_Box);
    Vec3D n2_1 = (pair.PT2.pl1)->GetMirrorLink()->GetTriangle()->GetNormalVector();
    Vec3D n2_2 = (pair.PT2.pl2)->GetMirrorLink()->GetTriangle()->GetNormalVector();
    Vec3D n2_3 = (pair.PT2.pl3)->GetMirrorLink()->GetTriangle()->GetNormalVector();
    if( n1.dot(n2,n2_1) < m_MinAngle ||
        n1.dot(n2,n2_2) < m_MinAngle ||
        n1.dot(n2,n2_3) < m_MinAngle ) {
        return false;
    }
    
    return true;
}
//--- v1,v2 and v3 may not have the correct trinagluation Orientation.
 //      >v         ;correct orientation means to disconnect links connected to vertices in the pot_triangle 2
//  l1  /  \>  l2   ;because for each triple of v1, v2 and v3, there is two ways to create a trinagle but each will
//    v1<---v3      ; leads to removal of diffierent edges.
bool Three_Edge_Scission::CorrectOrientation(pot_triangle &p1, pot_triangle &p2) {
    /**
     * @brief Ensures the correct triangulation orientation for the provided triangles.
     *
     * This function checks the orientation of the vertices in the `pot_triangle` p1 to ensure that
     * the v3 of each edge from p1 exists in the vertices of p2. If the orientation is incorrect,
     * it reverses it by swapping edges to mirror edges. If neither the original nor the mirrored links
     * have the correct orientation, it returns false.
     *
     * @param p1 The first potential triangle with three edges and three vertices.
     * @param p2 The second potential triangle with three vertices.
     * @return true if the orientation is correct or successfully corrected, false otherwise.
     */
    
    links* ml1 = (p1.pl1)->GetMirrorLink();
    if((p1.pl1)->GetV3() == p2.pv1 || (p1.pl1)->GetV3() == p2.pv2 || (p1.pl1)->GetV3() == p2.pv3){
        // Orientation is correct, no change is needed
        return true;
    }
    else if(ml1->GetV3() == p2.pv1 || ml1->GetV3() == p2.pv2 || ml1->GetV3() == p2.pv3){
        // Orientation is not correct, we reverse it
        links* ml2 = (p1.pl2)->GetMirrorLink();
        links* ml3 = (p1.pl3)->GetMirrorLink();
        vertex *v1 = p1.pv1;
        vertex *v2 = p1.pv2;
        p1.pv1 = v2;
        p1.pv2 = v1;
        
        p1.pl1 = ml1;
        p1.pl2 = ml3;
        p1.pl3 = ml2;


        return true;
    }
    else {  // the triangles are not connected
        return false;
    }
    
    return true;
}
void Three_Edge_Scission::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void Three_Edge_Scission::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
template<typename T>
bool Three_Edge_Scission::AddtoVectorCarefully(T* item, std::vector<T*>& vect) {
    // Check if the item already exists in the list
    for (typename std::vector<T*>::iterator it = vect.begin(); it != vect.end(); ++it) {
        if (*it == item)
            return false;
    }
    // Item does not exist, add it to the list
    vect.push_back(item);
    return true;
}
// this should be deleted at the end

template<typename T>
void Three_Edge_Scission::KeepOneOccurrence(std::vector<T*> &vec){
   
    std::sort(vec.begin(), vec.end()); // Sort the vector
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); // Remove duplicates
}
std::vector<links*> Three_Edge_Scission::GetEdgesWithInteractionChange(pair_pot_triangle &pair_t){
    
    /**
     * @brief Retrieve edges with interaction changes due to inclusions or vector fields.
     *
     * This function identifies and returns the edges (links) associated with the vertices of
     * two paired triangles that are expected to experience changes in interactions. The changes
     * in interactions could be due to inclusions owned by the vertices or the presence of vector
     * fields.
     *
     * The function operates as follows:
     * 1. It gathers all links associated with the vertices of the two input triangles (`pair_t`).
     * 2. It filters these links to include only those that either own inclusions or have associated
     *    vector fields.
     * 3. It removes duplicate links and links that are mirrors of each other.
     *
     * @param pair_t A reference to a pair of triangles, each represented by a `pair_pot_triangle` structure.
     *               The vertices of these triangles are used to gather and filter the links.
     * @return A vector of pointers to the `links` objects that have interaction changes.
     *
     * @note In development mode (`DEVELOPMENT_MODE == Enabled`),
     */
#if DEVELOPMENT_MODE == Enabled
    std::cout<<" DEVELOPMENT_MODE ID 9595473: This function can be  made  better\n ";
    // for example, precheck to select only links that both has includioon ...
#endif
    
    std::vector<links*> edge_with_interaction_change;
    std::vector<links *> all_temlinks;
//---> get all the links
    vertex *v11 = pair_t.PT1.pv1;
    vertex *v12 = pair_t.PT1.pv2;
    vertex *v13 = pair_t.PT1.pv3;
    vertex *v21 = pair_t.PT2.pv1;
    vertex *v22 = pair_t.PT2.pv2;
    vertex *v23 = pair_t.PT2.pv3;
    bool system_has_vf = false;
    if(v11->GetNumberOfVF() != 0){
        system_has_vf = true;
    }
    
    //=== inclusion interaction energy;
    if(v11->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v11->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
    }
    if(v12->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v12->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
    }
    if(v13->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v13->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
    }
    if(v21->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v21->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
    }
    if(v22->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v22->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
    }
    if(v23->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v23->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
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
std::string Three_Edge_Scission::CurrentState(){
    
    std::string state = AbstractDynamicTopology::GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+ Nfunction::D2S(m_Period);
    return state;
}


/*
 
 
 for (std::vector<triangle *>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it){
     (*it)->UpdateNormal_Area(&m_Box);   //trinagule normal and area should be obtained
 }
 for (std::vector<links *>::iterator it = m_pActiveL.begin() ; it != m_pActiveL.end(); ++it){
     (*it)->UpdateNormal();
     (*it)->UpdateShapeOperator(&m_Box);

 }
 for (std::vector<vertex *>::iterator it = m_pSurfV.begin() ; it != m_pSurfV.end(); ++it){
     (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(*it);
 }
 for (std::vector<links *>::iterator it = m_pActiveL.begin() ; it != m_pActiveL.end(); ++it){
     (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
 }
 */
