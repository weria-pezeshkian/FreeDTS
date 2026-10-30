#include <chrono>
#include <unordered_set>
#include <utility>
#include <ctime>
#include <set>
#include "TopologyChangeByTriangularPrism.h"
#include "State.h"
#include "MESH.h"
#include "vertex.h"
#include "./Registry/FactoryDynamicTopologyMethod.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */

TopologyChangeByTriangularPrism::TopologyChangeByTriangularPrism(std::string inputdata, State *pState) :
                m_pState(pState),
                m_StreamInputs (inputdata),
                m_PrismMapTopologyFile(""),
                m_pEdgeL(pState->GetMesh()->GetEdgeL()),
                m_pGhostL(pState->GetMesh()->GetGhostL()),
                m_pGhostT(pState->GetMesh()->GetGhostT()),
                m_pActiveT(pState->GetMesh()->GetActiveT()),
                m_pRightL(pState->GetMesh()->GetRightL()),
                m_pLeftL(pState->GetMesh()->GetLeftL()),
                m_pActiveL(pState->GetMesh()->GetActiveL()),
                m_pSurfV(pState->GetMesh()->GetSurfV()),
                m_pEdgeV(pState->GetMesh()->GetEdgeV()),
                m_Beta(pState->GetSimulation()->GetBeta()),
                m_DBeta(pState->GetSimulation()->GetDBeta()),
                m_MinLength2(pState->GetSimulation()->GetMinL2()),
                m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
                m_MinAngle(pState->GetSimulation()->GetMinAngle()),
                m_No_VectorFields_Per_V(pState->GetMesh()->GetNoVFPerVertex()),
                m_Box(pState->GetMesh()->Link2ReferenceBox()),
                m_pTriVoxelization(new Voxelization<triangle>())
{

}
TopologyChangeByTriangularPrism::~TopologyChangeByTriangularPrism(){
    delete m_pTriangularPrismBuilder ;
}
void TopologyChangeByTriangularPrism::Initialize() {


std::vector<std::string> input_data = Nfunction::Split(m_StreamInputs);
    if (input_data.empty()) {
        std::cerr << "---> Error: insufficient input data for '" << GetDefaultReadName() << "' command" << std::endl;
    }   
// Parse required parameter: Period
    m_Period = Nfunction::String_to_Int(input_data[0]);

// Parse optional parameter: PrismMapTopologyFile
    if (input_data.size() > 1 && !input_data[1].empty()) {
        m_PrismMapTopologyFile = input_data[1];
    } 
    else {
            std::cout << "---> Note: No topology prism file provided for '" 
              << GetDefaultReadName() << "' command. Using default behavior." << std::endl;
    }
    std::string note_txt = "This is a dynamic topology simulation with " + GetDefaultReadName() +" method.";
    Nfunction::ConsolePrint_Note(note_txt);
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;
    m_Surface_Genus = 1-(m_pSurfV.size()-m_pLeftL.size()+m_pActiveT.size())/2;
    
    m_pTriangularPrismBuilder = new TriangularPrismBuilder (&m_Box, m_MaxLength2, m_MinAngle, m_PrismMapTopologyFile);
}
//===========================================================
//========== MC MOVE: containing both fission and fusion move 
//==========================================================
bool TopologyChangeByTriangularPrism::MCMove(int step) {
    
    if(m_Period == 0 ||  step%m_Period != 0){
        return false;
    }
        
//========================================================
//--============== Fission
//======================================================
   std::vector<pair_pot_triangle> pair_list  = FindNecks();
    std::cout<<pair_list.size()<<" number of pot trinagles \n";
   /* if(pair_list.size() != 0) {// ScissionByMC
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(pair_list.size());
        pair_pot_triangle pair_T = pair_list[n];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(ScissionByMC(pair_T, thermal)){
            std::cout <<"accepted ScissionByMC "<<step<<"\n";
            m_AcceptedMoves++;
        }
    }*/ ///  if(pair_list.size() != 0) end ScissionByMC
    
//========================================================
//--============== should be removed
//======================================================
    /*
    std::clock_t start = std::clock();
     std::vector<fusion_site> all_possible_sites = FindPotentialFusionSites();
      std::clock_t end = std::clock();
     Calculate the duration
       double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        Print the time taken
       std::cout << "Time taken to execute FindPotentialFussionSites: " << duration << " seconds " << all_possible_sites.size() << std::endl;*/
//========================================================
//--============== Fusion
//======================================================         
  /*  std::vector<fusion_site> all_possible_sites = FindPotentialFusionSites();
    if(all_possible_sites.size() != 0) {// FussionByMove
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(all_possible_sites.size());
        fusion_site pair_T = all_possible_sites[n];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(FusionByMove(pair_T, thermal)){
            std::cout<<"accepted fusion "<<step<<"\n";
            std::cout <<"accepted fusion\n";
            m_AcceptedMoves++;
        }
    }*/ ///     if(all_possible_sites.size() != 0) {// end FussionByMove

    m_Surface_Genus = 1 - (m_pSurfV.size()-m_pLeftL.size()+m_pActiveT.size())/2;
    return true;
}
std::vector<fusion_site> TopologyChangeByTriangularPrism::FindPotentialFusionSites() {
    //////////////////////////////////////////////////////////////////////////////////////////
//  FindPotentialFusionSites
//
//  PURPOSE:
//  --------
//  Identify all valid potential fusion sites between triangles in the system.
//
//  OVERVIEW:
//  ---------
//  This function:
//   1. Builds a voxel grid representation of all active triangles.
//   2. Searches each triangle’s local 3x3x3 voxel neighborhood for candidates.
//   3. Filters candidate triangle pairs using:
//        - Periodic boundary condition distance cutoff
//        - Geometric exclusion rules (FusionSites_AreNotNeighbours)
//   4. Stores unique triangle pairs in a hash set to avoid duplicates.
//   5. Performs final validation on each pair to construct fusion_site objects.
//
//  RETURNS:
//  --------
//  A vector of validated fusion_site objects representing possible fusion events.
//
//////////////////////////////////////////////////////////////////////////////////////////
    
    
    std::unordered_set<std::pair<triangle*, triangle*>, PairHash> Possible_pairs;
    std::vector<fusion_site> Available_Sites;
    
    if (!VoxelizeTriangles(m_MaxLength2)) {
        std::cout << "Failed to make voxels\n";
        return Available_Sites;
    }

    for (auto it = m_pActiveT.begin(); it != m_pActiveT.end(); ++it)
    {
        triangle* t0 = *it;
        Vec3D T0center = t0->GetCentroid();

        Voxel<triangle>* pvoxel = t0->GetVoxel();
        if (!pvoxel) {
            std::cout << "Error -> NULL TRIANGLE VOXEL\n";
            std::abort();
        }

        for (int n = -1; n <= 1; n++)
        for (int m = -1; m <= 1; m++)
        for (int s = -1; s <= 1; s++)
        {
            Voxel<triangle>* pvox = pvoxel->GetANeighbourCell(n, m, s);
            if (!pvox) {
                std::cout << "Error -> NULL VOXEL DETECTED\n";
                std::abort();
            }

            const auto& VTri = pvox->GetContentObjects();

            for (auto itvox = VTri.begin(); itvox != VTri.end(); ++itvox)
            {
                triangle* t1 = t0;
                triangle* t2 = *itvox;

                if (t1 == t2) continue;

                Vec3D T1center = t2->GetCentroid();

                if (Nfunction::SquarePBCDistanceOfTwoPoint(
                        T0center, T1center, m_Box) >= m_MaxLength2)
                    continue;

                if (!FusionSites_AreNotNeighbours(t1, t2))
                    continue;

                if (t1 > t2) std::swap(t1, t2);  // to remove the repeated one

                Possible_pairs.insert({t1, t2});
            }
        }
    }

   // now we have unique pair of trinagles that are close
   // we can check them further
   //Performs final validation on each pair to construct fusion_site objects FusionSite_DistanceIsGood().
   for (const auto& p : Possible_pairs)
    {
        triangle* t1 = p.first;
        triangle* t2 = p.second;
        
        std::vector<TriangularPrism> PossiblePrism = m_pTriangularPrismBuilder->GeneratePossibleTopology(t1, t2);
        if(PossiblePrism.size()==0){
            continue;
        }
        for (const auto& topo : PossiblePrism){
            fusion_site p_T;
            p_T.t1 = t1;
            p_T.t2 = t2;
            p_T.topology  = topo;
            Available_Sites.push_back(p_T);
        }
    }
    
    return Available_Sites;
}
bool TopologyChangeByTriangularPrism::FusionSites_AreNotNeighbours(triangle *t1, triangle *t2){
 
 // this might be extended to avoid even second next neighbours 
    // ---- vertex sharing rejection ----
    vertex * v1 = t1->GetV1();
    vertex * v2 = t1->GetV2();
    vertex * v3 = t1->GetV3();

    vertex * u1 = t2->GetV1();
    vertex * u2 = t2->GetV2();
    vertex * u3 = t2->GetV3();

//-- check if they shapre any vertex. 
    if (v1 == u1 || v1 == u2 || v1 == u3 ||
        v2 == u1 || v2 == u2 || v2 == u3 ||
        v3 == u1 || v3 == u2 || v3 == u3)
    {
        return false;
    }
//-- check if they share any vertex with the vertex neighbour      
    std::vector <vertex *> pNv1 = v1->GetVNeighbourVertex();
    for (auto it : pNv1) {
        if (it == u1 || it == u2 || it == u3){
            return false;
        }
    }
    
    std::vector <vertex *> pNv2 = v2->GetVNeighbourVertex();
    for (auto it : pNv2) {
        if (it == u1 || it == u2 || it == u3){
            return false;
        }
    }
    
    std::vector <vertex *> pNv3 = v3->GetVNeighbourVertex();
    for (auto it : pNv3) {
        if (it == u1 || it == u2 || it == u3){
            return false;
        }
    }
    
    
    return true;
        
}
//========================================================================
//=====================  Fusion function =================================
//========================================================================
bool TopologyChangeByTriangularPrism::FusionByMove(fusion_site &pair_tri, double thermal){
    
    // There is a lot to be added to this function
       // if(m_AcceptedMoves!=0)
         //return false;

    double new_energy = 0;
    double old_energy = 0;
    
    triangle* t1 = pair_tri.t1;
    triangle* t2 = pair_tri.t2;
    
    vertex * v1 = t1->GetV1();
    vertex * v2 = t1->GetV2();
    vertex * v3 = t1->GetV3();

    vertex * u1 = t2->GetV1();
    vertex * u2 = t2->GetV2();
    vertex * u3 = t2->GetV3();
    
    std::vector<vertex*> Vver;
    Vver.insert(Vver.end(),{v1, v2, v3, u1, u2, u3});
    
    // get old energy
    for (const auto& pV : Vver){
        old_energy += pV->GetEnergy();
       // old_energy += pV->GetBindingEnergy(); // vector field
    }
    
    //-- get the energy for interaction
        std::vector<links*> Affected_links  = Get_EdgesFusionAffect(Vver);

        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
                (*it)->Copy_InteractionEnergy();
                (*it)->Copy_VFInteractionEnergy();
                old_energy += 2 * (*it)->GetIntEnergy();
                old_energy += 2 * (*it)->GetVFIntEnergy();
        }
        
        // and more terms: global variables 
            // Obtain and sum the initial global variables that might change
            double old_Tvolume = 0.0, old_Tarea = 0.0, old_Tcurvature = 0.0;
            double new_Tvolume = 0.0, new_Tarea = 0.0, new_Tcurvature = 0.0;
            /*if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
                m_pState->GetVAHGlobalMeshProperties()->CalculateALinkTrianglesContributionToGlobalVariables(p_edge, old_Tvolume, old_Tarea, old_Tcurvature);
            }*/
        //====
           // we know perform fusion and our mesh is no longer the old one.
           // this will also updates trinagule and links normal and shape operator 

             fusion_outcome fusion_data;          
            if (!Fuse_MeshViaTwoTriangles(pair_tri, fusion_data)) {
                std::cerr << "ERROR: Fuse_MeshViaTwoTriangles() returned false. "
                        "Fusion was expected to succeed for the selected triangle pair."
                << std::endl;
            }
        
            // update curvature of the vertices
            for (const auto& pV : Vver){
                (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(pV);

            }
            for (const auto& pV : Vver){
                new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(pV);
                //new_energy += pV->GetBindingEnergy(); // vector field
            }
            

            for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
                    (*it)->Copy_InteractionEnergy();
                    (*it)->Copy_VFInteractionEnergy();
                    new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
                    
                if(v1->GetNumberOfVF() != 0 ){
                   for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                     new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                    }    
                }
            }
            // new links need to be updated
            for (std::vector<links *>::iterator it = (fusion_data.pHalfnewLinks).begin() ; it != (fusion_data.pHalfnewLinks).end(); ++it){
                    new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
                    
                if(v1->GetNumberOfVF() != 0 ){
                   for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                     new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                    }    
                }
            }
            
        

    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy ;



   /* double energy0 = m_pState->GetEnergyCalculator()->GetEnergy();
    double Final_energy = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();    
    std::cout<<"f energy: "<<Final_energy<<" 0energy "<<   energy0+diff_energy <<"\n";
    */
    double U = m_Beta * tot_diff_energy - m_DBeta;
    //---> accept or reject the move
    if(2>1){//U <= 0 || exp(-U) > thermal ) {
        m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
        //std::cout<<"accepted \n";
        return true;
    }
    else {
        
        if(!Reverse_Fuse_MeshViaTwoTriangles(fusion_data)){
                Nfunction::ConsolePrint_Error("-- error should not happen ");
        }
        
        for (auto* li : Affected_links) {
            li->Reverse_InteractionEnergy();
            li->Reverse_VFInteractionEnergy();
        }
        for (auto* li : fusion_data.newLinks ) {
            li->Reverse_InteractionEnergy();
            li->Reverse_VFInteractionEnergy();
        }
        return false;
    }
        
        
        
        
    
    return false;
}

std::vector<links*> TopologyChangeByTriangularPrism::Get_EdgesFusionAffect(std::vector<vertex*> &Vver){  
/**
 * @brief Collects unique undirected edges from a set of vertices
 * 
 * Iterates through all links incident to the given vertices and returns
 * each edge only once, treating (v1→v2) and (v2→v1) as the same edge.
 * Uses hash set for O(N) average time complexity.
 * 
 * @param Vver Vector of vertex pointers to collect edges from
 * @return Vector of unique link pointers (one per undirected edge)
 * 
 * @note Order of returned links is not guaranteed
 * @note Result contains no duplicates even if two vertices share the same edge
 */

  
    auto hash = [](const std::pair<vertex*, vertex*>& p) {
        return std::hash<vertex*>()(p.first) ^ (std::hash<vertex*>()(p.second) << 1);
    };
    
    std::unordered_set<std::pair<vertex*, vertex*>, decltype(hash)> seen_edges(0, hash);
    std::vector<links*> result;
    result.reserve(100);
    seen_edges.reserve(100);
    
    for (auto* pv : Vver) {
        for (auto* link : pv->GetVLinkList()) {
            auto v1 = link->GetV1();
            auto v2 = link->GetV2();
            auto key = std::make_pair(std::min(v1, v2), std::max(v1, v2));
            
            if (seen_edges.insert(key).second) {
                result.push_back(link);
            }
        }
    }
    return result;
}
//========================================================================
//=====================  Fusion function =================================
//========================================================================
bool TopologyChangeByTriangularPrism::Reverse_Fuse_MeshViaTwoTriangles(fusion_outcome &fusion_mesh){
    /**
 * @brief Reverses a previous mesh fusion operation, restoring the original topology
 * 
 * This function is the inverse of Fuse_MeshViaTwoTrinagles(). It restores the mesh
 * to its state before the fusion by moving triangles and links back from active
 * to ghost pools and reverting all modified objects to their saved copies.
 * 
 * @param fusion_mesh   A fusion_outcome structure containing all elements created
 *                      during the forward fusion operation, including:
 *                      - t1, t2: The original two triangles that were fused
 *                      - newTriangles: The 6 triangles created during fusion
 *                      - newLinks: The 12 links created during fusion
 *                      - oldLinks: The original links that were modified
 *                      - vertices: The 6 vertices involved in the fusion
 * 
 * @return true         Always returns true on successful reversal
 * 
 * @details The reversal process consists of six steps:
 * 
 * 1. **Ghost Triangle Cleanup**: Removes the two original triangles (t1, t2)
 *    from the front of the ghost pool where they were placed during fusion.
 * 
 * 2. **Active Triangle Restoration**: 
 *    - Removes the 6 newly created triangles from the active pool
 *    - Restores the original two triangles (t1, t2) back to the active pool
 *    - Moves the 6 new triangles to the ghost pool for future use
 * 
 * 3. **Link Pool Restoration**:
 *    - Removes the 12 new links from the active link pool
 *    - Removes the 6 pairs from the left/right link containers
 *    - Moves the 12 new links to the ghost link pool
 * 
 * 4. **Vertex State Restoration**: Calls Reverse2PreviousCopy() on all 6 vertices
 *    to restore their saved state (triangle lists, coordinates, etc.)
 * 
 * 5. **Link State Restoration**: Calls Reverse2PreviousCopy() on all modified
 *    original links to restore their neighbor relationships and geometric properties
 * 
 * 6. **Ghost Pool Management**: Appends the 12 new links to the ghost link pool
 *    for future reuse in subsequent fusion operations
 * 
 * @warning This function assumes that fusion_mesh contains valid pointers and that
 *          the ghost pools have sufficient capacity. No error checking is performed.
 * 
 * @note The order of operations is critical. Pools are resized before pushing
 *       back elements to maintain proper container sizes and avoid memory issues.
 * 
 * @see Fuse_MeshViaTwoTrinagles() The forward fusion operation that this reverses
 */
    
    // remove the first two members. Because we added t1 and t2 to the first member 
    if (m_pGhostT.size() >= 2) {
        m_pGhostT.erase(m_pGhostT.begin(), m_pGhostT.begin() + 2);
    }
    // remove the new trinagle from active
    m_pActiveT.resize(m_pActiveT.size() - 6);
    m_pActiveT.push_back(fusion_mesh.t1);
    m_pActiveT.push_back(fusion_mesh.t2);
    // add back the new trinagle to the ghost
    for(auto* tri: fusion_mesh.newTriangles){
        m_pGhostT.push_back(tri);
    }
    m_pActiveL.resize(m_pActiveL.size() - 12);
    m_pRightL.resize(m_pRightL.size() - 6);
    m_pLeftL.resize(m_pLeftL.size() - 6);
    
    
    for(auto* pver: fusion_mesh.vertices){
        pver->Reverse2PreviousCopy();
    }
    for(auto* li: fusion_mesh.oldLinks){
        li->Reverse2PreviousCopy();
    }
    for(auto* li: fusion_mesh.newLinks){
        m_pGhostL.push_back(li);
    }

    
 return true;   
}
//========================================================================
//=====================  Fusion function =================================
//========================================================================
bool TopologyChangeByTriangularPrism::Fuse_MeshViaTwoTriangles(fusion_site &pair_tri, fusion_outcome &outcome){
// -----------------------------------------------------------------------------
// bool TopologyChangeByTriangularPrism::Fuse_MeshViaTwoTrinagles(fusion_site &pair_tri)
// -----------------------------------------------------------------------------
// Fuses two triangles into a new prism-like local topology using ghost
// triangle/link pools. The function:
//   1. Validates available ghost resources
//   2. Removes old triangles from active mesh
//   3. Builds 6 new triangles from topology mapping
//   4. Creates and connects new links
//   5. Assigns mirror relationships between links
//   6. Recomputes normals and shape operators
// -----------------------------------------------------------------------------

         
if (m_pGhostT.size() < 6 || m_pGhostL.size() < 12) {
    std::cout << "---> Warning: Ghost counters are full. Fusion will not be performed. "
              << "This is a code limitation, not a physics limitation. "
              << "Please restart the simulation to restore correct physics behavior.\n";
    return false;
}
    
    triangle* t1 = pair_tri.t1;
    triangle* t2 = pair_tri.t2;
    RemoveFromTriangleList(t1, m_pActiveT);
    RemoveFromTriangleList(t2, m_pActiveT);
    // we add them to the begining of the ghost, so when we take out from the end, for the moment, we do not affect these objects
    m_pGhostT.insert(m_pGhostT.begin(), t1);
    m_pGhostT.insert(m_pGhostT.begin(), t2);

        
    vertex * v1 = t1->GetV1();
    vertex * v2 = t1->GetV2();
    vertex * v3 = t1->GetV3();
    vertex * u1 = t2->GetV1();
    vertex * u2 = t2->GetV2();
    vertex * u3 = t2->GetV3();
    std::vector<vertex*> Vver;
    Vver.insert(Vver.end(),{v1, v2, v3, u1, u2, u3});
    for (auto* ver: Vver){
        
        if(!ver->SetCopy()){
            Nfunction::ConsolePrint_Error(" Error-> copying faild ");
        }
    }
    v1->RemoveFromTraingleList(t1);
    v2->RemoveFromTraingleList(t1);
    v3->RemoveFromTraingleList(t1);
    u1->RemoveFromTraingleList(t2);
    u2->RemoveFromTraingleList(t2);
    u3->RemoveFromTraingleList(t2);
    

    std::vector<triple> VTriples = (pair_tri.topology).VTriples;
    std::vector<triangle*> pnewTrinagles;
    std::vector<links*> pnewLinks;
    std::vector<links*> pHalfnewLinks;
    std::vector<links*> poldLinks;

    for (int i= 0; i < 6; i++){
            triangle* prism_t = m_pGhostT.back();
            m_pGhostT.pop_back();
            m_pActiveT.push_back(prism_t);
            const int id1 = VTriples[i][0];
            const int id2 = VTriples[i][1];
            const int id3 = VTriples[i][2];
            vertex * tv1 = Vver[id1];
            vertex * tv2 = Vver[id2];
            vertex * tv3 = Vver[id3];

            prism_t->UpdateVertex(tv1, tv2, tv3);
            pnewTrinagles.push_back(prism_t);
            tv1->AddtoTraingleList(prism_t);
            tv2->AddtoTraingleList(prism_t);
            tv3->AddtoTraingleList(prism_t);
            
            // lets now create the links 
            links* old_link = tv1->GetConnectingLink(tv2);
            if(!old_link->SetCopy()){
                Nfunction::ConsolePrint_Error(" Error-> copying link faild ");
            }
            links* prism_l1 = m_pGhostL.back();
            m_pGhostL.pop_back();
            m_pActiveL.push_back(prism_l1);
            links* prism_l2 = m_pGhostL.back();
            m_pGhostL.pop_back();
            m_pActiveL.push_back(prism_l2);
            
            prism_l1->UpdateV(tv2, tv3, tv1);
            prism_l2->UpdateV(tv3, tv1, tv2);
            prism_l1->UpdateMirrorFlag(false);  //     // since the  mirrors have not been found yet
            prism_l2->UpdateMirrorFlag(false);

            prism_l1->UpdateNeighborLink1(prism_l2);
            prism_l1->UpdateNeighborLink2(old_link);
            prism_l2->UpdateNeighborLink1(old_link);
            prism_l2->UpdateNeighborLink2(prism_l1);
            old_link->UpdateNeighborLink1(prism_l1);
            old_link->UpdateNeighborLink2(prism_l2);
            
            tv2->AddtoLinkList(prism_l1);
            tv3->AddtoLinkList(prism_l2);
            
            old_link->UpdateTriangle(prism_t);
            prism_l1->UpdateTriangle(prism_t);
            prism_l2->UpdateTriangle(prism_t);
            
            poldLinks.push_back(old_link);
            pnewLinks.push_back(prism_l1);
            pnewLinks.push_back(prism_l2);

    }

    for (size_t i = 0; i < pnewLinks.size(); i++) {
            
            links* l1 = pnewLinks[i];
        if (l1->GetMirrorFlag()) {
            continue;
        }

        links* l2 = nullptr;
        bool found = false;

        for (size_t j = i + 1; j < pnewLinks.size(); j++) {
            l2 = pnewLinks[j];
            if (l2->GetMirrorFlag()) {
                continue;
            }
            if (l1->GetV1() == l2->GetV2() && l1->GetV2() == l2->GetV1()){
                found = true;
                break;
            }
        }

        if (!found) {
              std::cerr << "\033[1;31m"  // Bold red
              << "╔══════════════════════════════════════════════════════════════╗\n"
              << "║                    !!! FATAL ERROR !!!                       ║\n"
              << "╚══════════════════════════════════════════════════════════════╝\n"
              << "\033[0m"  // Reset
              << "\033[1;31m"   // Blue
              << "---> In "<<TopologyChangeByTriangularPrism::GetDefaultReadName()<<" command, "
              << " function Fuse_MeshViaTwoTrinagles(), \n"
              << "      we are matching mirror edges, but one found that does not have any  \n"
              << "      links vertices id  "<<pnewLinks[i]->GetV1()->GetVID()<<"   "<<pnewLinks[i]->GetV2()->GetVID()<<"\n"

              << "══════════════════════════════════════════════════════════════\n"
              << "\033[0m";  // Reset
              return false;
             }
        m_pRightL.push_back(l1);
        m_pLeftL.push_back(l2);
        pHalfnewLinks.push_back(l1);
        l1->UpdateMirrorLink(l2);
        l2->UpdateMirrorLink(l1);
        l1->UpdateMirrorFlag(true);
        l2->UpdateMirrorFlag(true);
    }

    // add the vertices to the nighbour list
    for (auto li : pHalfnewLinks){
        li->GetV1()->AddtoNeighbourVertex(li->GetV2());
        li->GetV2()->AddtoNeighbourVertex(li->GetV1());
    }

   // update geometry of the trinagles and the links
    for (auto tri : pnewTrinagles){
        tri->UpdateNormal_Area(&m_Box);
    }
    for (auto le : pHalfnewLinks){
        le->UpdateNormal();
        le->UpdateShapeOperator(&m_Box);
    }
    for (auto le : poldLinks){
        le->UpdateNormal();
        le->UpdateShapeOperator(&m_Box);
    }
/*int nn =0;
for (auto ver : Vver) {
    
   inclusion* inc = ver->GetInclusion();
   if(nn<3){
   inc->m_IncType = m_pState->GetMesh()->GetInclusionType()[2];
   }
   else
   {
          inc->m_IncType = m_pState->GetMesh()->GetInclusionType()[3];

   }
nn++;
}
*/

    outcome.newTriangles = pnewTrinagles;
    outcome.newLinks = pnewLinks;
    outcome.oldLinks = poldLinks;
    outcome.pHalfnewLinks = pHalfnewLinks;
    outcome.vertices = Vver;
    outcome.t1 = t1;
    outcome.t2 = t2;

    return true;
}

bool TopologyChangeByTriangularPrism::ScissionByMC(pair_pot_triangle &pair_t, double thermal){
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
    
    return false;
}
// it creates a triangle and place it to the active trinagles list
triangle * TopologyChangeByTriangularPrism::CreateATriangleFromAPotentialTriangle(pot_triangle &p1) {
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
std::vector <triangle *> TopologyChangeByTriangularPrism::DoAScission(pair_pot_triangle &pair){
    
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
bool TopologyChangeByTriangularPrism::ReverseAScission(pair_pot_triangle &pair , triangle *t1, triangle *t2){

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
void TopologyChangeByTriangularPrism::MakeALinkGhost(links *p_links){
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
void TopologyChangeByTriangularPrism::MakeATriangleGhost(triangle *p_tri){

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
void TopologyChangeByTriangularPrism::ActivateAGhostTriangle(triangle *p_tri){

    RemoveFromTriangleList(p_tri,m_pGhostT);
    AddtoVectorCarefully(p_tri, m_pActiveT);
    //--- add the triangles to the vertex list
    p_tri->GetV1()->AddtoTriangleListCarefully(p_tri);
    p_tri->GetV2()->AddtoTriangleListCarefully(p_tri);
    p_tri->GetV3()->AddtoTriangleListCarefully(p_tri);

    return;
}
void TopologyChangeByTriangularPrism::ActivateAGhostLink(links *p_links){
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
auto canonicalize = [](vertex*& a, vertex*& b, vertex*& c) {
    if (a > b) std::swap(a, b);
    if (a > c) std::swap(a, c);
    if (b > c) std::swap(b, c);
};
//== this function get a mesh and search through it and finds possible fission sites
std::vector<pair_pot_triangle> TopologyChangeByTriangularPrism::FindNecks() {
    /**
 * @brief Finds unique potential triangles ("triangle frames") in the surface mesh and groups them into neck candidates.
 *
 * This function scans the surface vertex connectivity to identify all valid triplets of vertices
 * (v1, v2, v3) that are pairwise connected by edges but do not yet form an explicit triangle in the mesh.
 *
 * To avoid redundant entries, each triangle is first canonicalized (sorted by vertex pointer address)
 * and then inserted into a hash set to ensure uniqueness regardless of traversal order or permutation.
 *
 * After constructing the set of unique potential triangles, the function checks all pairs of these
 * triangles to determine whether they form a valid "neck" configuration using Is_A_Neck().
 *
 * @return A vector of paired potential triangles representing detected neck structures in the mesh.
 */
    std::vector<pot_triangle> trinagle_frame;
    std::vector<pair_pot_triangle> pair_list;

    std::unordered_set<PotTriangleKey, PotTriangleKeyHash> unique_tris;

    int id = 0;

    for (auto pv1 : m_pSurfV) {

        auto nv = pv1->GetVNeighbourVertex();

        for (size_t i = 0; i < nv.size(); i++) {
            for (size_t j = i + 1; j < nv.size(); j++) {

                vertex* pv2 = nv[i];
                vertex* pv3 = nv[j];

                if (!pv1->IsAnyTriangle(pv2, pv3)) {  // vertices are nighbour of pv1 but not forming a triangle already
                    if(!pv2->IsThereAConnectingLink(pv3)) // pv2-pv3 are connected
                            continue;
                    // ---- canonical ordering
                    canonicalize(pv1, pv2, pv3);

                    // ---- uniqueness check
                    PotTriangleKey key{pv1, pv2, pv3};
                    if (unique_tris.insert(key).second) {

                        pot_triangle PotT;
                        PotT.id = id++;
                        PotT.pv1 = pv1;
                        PotT.pv2 = pv2;
                        PotT.pv3 = pv3;

                        trinagle_frame.push_back(PotT);
                    }
                }
            }
        }
    }

    int idpair = 0;
    for (int i = 0; i < (int)trinagle_frame.size(); i++) {
        for (int j = i + 1; j < (int)trinagle_frame.size(); j++) {

            pair_pot_triangle temp;

            if (Is_A_Neck(trinagle_frame[i], trinagle_frame[j], temp)) {
                temp.id = idpair++;
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
bool TopologyChangeByTriangularPrism::Is_A_Neck(pot_triangle potT1, pot_triangle potT2, pair_pot_triangle & neck) {
    /**
     * @brief Checks if two potential triangles (pot_triangle) are well-connected to form a neck for fission.
*/

neck.state = false;
    
vertex* v_ver[3] = {potT1.pv1, potT1.pv2, potT1.pv3};
vertex* u_ver[3] = {potT2.pv1, potT2.pv2, potT2.pv3};

int i = 0, j = 0;

//---> if the two potential trinagles share a vertex, then they should be discarded
    // Ensure the two potential triangles do not share any vertices.
while (i < 3 && j < 3) {
    if (v_ver[i] == u_ver[j]) {
        return false;
    }
    (v_ver[i] < u_ver[j]) ? ++i : ++j;
}
    
 //--    
 //----> this check if the potential triangles are connected or not
 std::vector<links*> bride_links;  // we collect all the links that can be removed if the neck is cut., we have to remove their mirror too
int connections[3][3] = {};
    for (int i = 0; i<3; i++){
        for (int j=0; j<3; j++){
            links* clink;
            if(v_ver[i]->IsThereAConnectingLink(u_ver[j], clink)){
                connections[i][j] = 0;
                bride_links.push_back(clink);
            }
        }  
    }
// Check if any row or column has sum of zero; this means it is well connected as a neck
for (int i = 0; i < 3; i++) {
    int rowSum = 0, colSum = 0;
    for (int j = 0; j < 3; j++) {
        rowSum += connections[i][j];
        colSum += connections[j][i];
    }
    if (rowSum == 0 || colSum == 0) return false;
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
bool TopologyChangeByTriangularPrism::CheckFaceAngle(pair_pot_triangle &pair){
    
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
bool TopologyChangeByTriangularPrism::CorrectOrientation(pot_triangle &p1, pot_triangle &p2) {
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
void TopologyChangeByTriangularPrism::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void TopologyChangeByTriangularPrism::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void TopologyChangeByTriangularPrism::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void TopologyChangeByTriangularPrism::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
template<typename T>
bool TopologyChangeByTriangularPrism::AddtoVectorCarefully(T* item, std::vector<T*>& vect) {
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
void TopologyChangeByTriangularPrism::KeepOneOccurrence(std::vector<T*> &vec){
   
    std::sort(vec.begin(), vec.end()); // Sort the vector
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); // Remove duplicates
}
std::vector<links*> TopologyChangeByTriangularPrism::GetEdgesWithInteractionChange(pair_pot_triangle &pair_t){
    
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
bool TopologyChangeByTriangularPrism::VoxelizeTriangles(double voxelsize) {
// Builds a voxelization of all active triangles using their PBC-correct centroids.
// Centroids are updated before insertion to ensure consistent spatial hashing.
// Returns a fully constructed voxel grid for triangle lookup / interaction queries.


    m_pTriVoxelization->SetBox(&m_Box);
    m_pTriVoxelization->UpdateVoxelSize(voxelsize, voxelsize, voxelsize);

    // Only keep this if it is actually needed separately from UpdateVoxelSize
    m_pTriVoxelization->SetDefaultVoxelSize(voxelsize, voxelsize, voxelsize);

    // Recompute centroids (PBC-correct)
    for (triangle* t : m_pActiveT)
    {
        t->CalculateCentroid(m_Box);
    }

    m_pTriVoxelization->Voxelize(m_pActiveT);

    return true;
}
std::string TopologyChangeByTriangularPrism::CurrentState(){
    
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
/*
    Static registration for TopologyChangeByTriangularPrism
    Automatically registers the class with FactoryDynamicTopologyMethod
    so it can be created dynamically from input streams.
*/
//
class RegistryTopologyChangeByTriangularPrism
{
public:

    RegistryTopologyChangeByTriangularPrism()
    {
        FactoryDynamicTopologyMethod::Instance().Register(
            TopologyChangeByTriangularPrism::GetDefaultReadName(),
            Create);
    }

private:

    static AbstractDynamicTopology* Create(
        std::istream& input,
        State* state)
    {
        //int period = 0;
       // input >>  period;
        
        std::string inputdata;
        std::getline(input, inputdata);   // consume rest of line
        
        return new TopologyChangeByTriangularPrism(
            inputdata,
            state);   // state is always last
    }
};
/*
---------------------------------------------------------------
Static Registration Object
---------------------------------------------------------------
*/
namespace
{
    RegistryTopologyChangeByTriangularPrism
        register_TopologyChangeByTriangularPrism;
}
