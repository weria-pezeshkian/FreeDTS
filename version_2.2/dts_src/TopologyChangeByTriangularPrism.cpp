/**
 * @class TopologyChangeByTriangularPrism
 * @brief Monte Carlo dynamic topology algorithm for membrane fusion and fission
 *        based on local triangular-prism transformations.
 *
 * This class implements topology-changing Monte Carlo moves for triangulated
 * membrane surfaces. Fusion and fission events are represented as local
 * transformations of a six-vertex triangular-prism structure, allowing the
 * membrane connectivity and topology to evolve while maintaining a valid
 * oriented triangulation.
 *
 * Candidate topology changes are generated dynamically:
 * - Potential fission sites are identified from neck-like structures.
 * - Potential fusion sites are identified from spatially separated nearby
 *   triangle pairs.
 *
 * Each proposed move is accepted or rejected according to the Metropolis-
 * Hastings criterion. The proposal probability accounts for the number of
 * available reverse moves, ensuring detailed balance when the number of
 * possible fusion and fission sites changes after a successful topology
 * transformation.
 *
 * The Gaussian curvature energy is defined as:
 *
 *        E_G = K_G \int K dA
 *
 * where K_G follows the convention used in the literature. Since the Gaussian
 * curvature integral is a topological quantity for closed surfaces, fusion
 * and fission moves include the corresponding topological energy contribution.
 *
 * The total energy difference used in the acceptance criterion may include:
 * - Local membrane bending energy
 * - Inclusion interaction energy
 * - Vector-field interaction energy
 * - Volume constraint energy
 * - Total area constraint energy
 * - Global curvature energy
 * - Gaussian curvature energy
 *
 * The algorithm supports:
 * - Detection of admissible fusion and fission configurations.
 * - Local topology modification using triangular-prism transformations.
 * - Exact reversal of rejected Monte Carlo moves.
 * - Updating of mesh connectivity, curvature, and interaction energies.
 * - Dynamic calculation of surface genus from the Euler characteristic.
 *
 * Memory management is optimized through the use of preallocated ghost
 * triangles and ghost links, avoiding dynamic allocation during topology
 * changes. Rejected moves are restored using stored local copies of affected
 * vertices, links, and interaction states.
 *
 * Assumptions:
 * - The membrane is represented as an oriented triangular mesh.
 * - Local topology transformations preserve mesh validity and orientation.
 * - Sufficient ghost elements are available for proposed transformations.
 *
 * @see TriangularPrismBuilder
 * @see AbstractDynamicTopology
 */

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
                m_pTriVoxelization(new Voxelization<triangle>()),
                m_No_Vectorfield(0),
                m_Local_KG(0)
{

}
TopologyChangeByTriangularPrism::~TopologyChangeByTriangularPrism(){
    delete m_pTriangularPrismBuilder ;
}
void TopologyChangeByTriangularPrism::Initialize() {


std::vector<std::string> input_data = Nfunction::Split(m_StreamInputs);
    if (input_data.size() < 2 ) {
        std::cerr << "---> Error: insufficient input data for '" << GetDefaultReadName() << "' command" << std::endl;
    }   
// Parse required parameter: Period
    m_Period = Nfunction::String_to_Int(input_data[0]);
    m_Local_KG = Nfunction::String_to_Double(input_data[1]);
    if (m_Period < 0) {
        throw std::invalid_argument(
            "Period must be greater than zero for '" +
            GetDefaultReadName() +
            "' command."
        );
    }

    if (m_Local_KG < 0.0) {
        Nfunction::ConsolePrint_Warning(
            "The input k_G is negative. Note that this command uses "
            "E_b = k_G * integral(K dA)."
        );
    }
// Parse optional parameter: PrismMapTopologyFile
    if (input_data.size() > 2 && !input_data[2].empty()) {
        m_PrismMapTopologyFile = input_data[2];
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

    // number of vf in the surface
    m_No_Vectorfield = m_pSurfV[0]->GetNumberOfVF();

//========================================================
//--============== Fission
//======================================================
    std::vector<fission_site> All_fission_Sites  = FindNecks();
    if(All_fission_Sites.size() != 0) {// ScissionByMC
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(All_fission_Sites.size());
        fission_site f_site = All_fission_Sites[n];    
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(ScissionByMC(f_site, thermal, All_fission_Sites.size())){
            m_AcceptedMoves++;
        }
    } ///  if(pair_list.size() != 0) end ScissionByMC
    
//========================================================
//--============== Fusion
//======================================================         
    std::vector<fusion_site> all_possible_sites = FindPotentialFusionSites();
    if(all_possible_sites.size() != 0) {// FussionByMove
        int n = m_pState->GetRandomNumberGenerator()->IntRNG(all_possible_sites.size());
        fusion_site pair_T = all_possible_sites[n];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        m_NumberOfAttemptedMoves++;
        if(FusionByMove(pair_T, thermal, all_possible_sites.size())){
            m_AcceptedMoves++;
        }
    } ///     if(all_possible_sites.size() != 0) {// end FussionByMove

    m_Surface_Genus = 1 - (m_pSurfV.size()-m_pLeftL.size()+m_pActiveT.size())/2;
    return true;
}
bool TopologyChangeByTriangularPrism::ScissionByMC(fission_site &f_site, double thermal, int n_fi_sites){
/**
 * @brief Attempts a Monte Carlo scission move by splitting a triangular-prism neck region.
 *
 * This routine proposes a topological scission at the specified fission site,
 * evaluates the resulting energy change, and accepts or rejects the move using
 * the Metropolis criterion.
 *
 * The algorithm:
 * - Stores the current local state and interaction energies.
 * - Computes the energy contribution of all vertices, links, and global
 *   geometric terms affected by the scission.
 * - Performs the trial scission, replacing the prism neck with two separated
 *   triangular caps and updating local connectivity.
 * - Recomputes the affected energies and global constraint contributions.
 * - Evaluates the total energy difference and applies the Metropolis test.
 * - If rejected, restores all topology, interaction energies, and copied state.
 *
 * Global energy terms that may contribute to the acceptance probability include:
 * - Volume constraint energy
 * - Total area constraint energy
 * - Global curvature energy
 *
 * @param f_site
 *        Description of the candidate scission region, including the prism
 *        vertices, links, and neighboring topology required to perform and
 *        potentially reverse the move.
 *
 * @param thermal
 *        Uniform random number in the range [0,1) used in the Metropolis
 *        acceptance test.
 *
 * @return true if the scission move is accepted and the new topology is kept;
 *         false if the move is rejected and the original configuration is
 *         restored.
 */
    // Check if there are enough links and triangles in the repository
    if (m_pGhostT.size() < 4 || m_pGhostL.size() < 4) {
        std::string message = " --->note: the number of the links and triangles in the repository is not enough, restart the simulations \n";
        *(m_pState->GetTimeSeriesLog())<<message;
        std::cout << message ;
        return false;
    }
    
    double new_energy = 0;
    double old_energy = 0;

    
//---> calculate old energies    
//---- effected vertex energy
    for (auto it_v : f_site.v_ver){
        it_v->SetCopy();
        old_energy += it_v->GetEnergy();
        old_energy += it_v->GetBindingEnergy(); // vector field

    }
    for (auto it_v : f_site.u_ver){
        it_v->SetCopy();
        old_energy += it_v->GetEnergy();
        old_energy += it_v->GetBindingEnergy(); // vector field
    }
    for (auto it_l : f_site.V_links){
        it_l->SetCopy();
         it_l->Copy_VFInteractionEnergy();
    }
    for (auto it_l : f_site.U_links){
        it_l->SetCopy();
         it_l->Copy_VFInteractionEnergy();
    }
//---- edge effects    
    std::vector<links*> Affected_links_old = GetEdgesWithInteractionChange(f_site);
    for (std::vector<links *>::iterator it = Affected_links_old.begin() ; it != Affected_links_old.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }    
            // and more terms: global variables 
            // Obtain and sum the initial global variables that might change
            double old_Tvolume = 0.0, old_Tarea = 0.0, old_Tcurvature = 0.0;
            double new_Tvolume = 0.0, new_Tarea = 0.0, new_Tcurvature = 0.0;
            if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
                std::vector<vertex*> verts;
                verts.insert(verts.end(), f_site.v_ver.begin(), f_site.v_ver.end());
                verts.insert(verts.end(), f_site.u_ver.begin(), f_site.u_ver.end());
                m_pState->GetVAHGlobalMeshProperties()->CalculateAPrismFacesContributionToGlobalVariables(verts, f_site.C_Links, old_Tvolume, old_Tarea, old_Tcurvature);
            }
            
//======== new mesh            
                    std::vector<triangle *> new_triangles = DoAScission(f_site);
//======== new mesh 

    for (auto it_v : f_site.v_ver){
        new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(it_v);
    }
    for (auto it_v : f_site.u_ver){
        new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(it_v);
    }
    
    std::vector<links*> Affected_links_new = GetEdgesWithInteractionChange(f_site);
    for (std::vector<links *>::iterator it = Affected_links_new.begin() ; it != Affected_links_new.end(); ++it){
            new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);

            if(m_No_Vectorfield != 0 ){
                for( int vf_layer = 0; vf_layer< m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++){
                new_energy +=  (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(vf_layer, *it);
                }
            }
    }
    
            // and more terms: global variables 
            // Obtain and sum the initial global variables that might change
            if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
                m_pState->GetVAHGlobalMeshProperties()->CalculateAPrismBasesContributionToGlobalVariables(new_triangles[0], new_triangles[1], new_Tvolume, new_Tarea, new_Tcurvature);
            }
     
        //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
    

    //==== MC
    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy + dE_volume + dE_t_area + dE_g_curv;
    
    // gaussian energy contribution
    tot_diff_energy += m_Local_KG * 4 * PI;

    
    
    double U = m_Beta * tot_diff_energy - m_DBeta;
    
    
    
    std::vector<fusion_site> fu_sites = FindPotentialFusionSites(); // this gives us N+S
    int n_fu_sites = fu_sites.size();
    // n_fi_sites; is the number of fission site before the move (M)
    if(n_fu_sites == 0){
        std::cout<<"---> this should not happen 2023342 \n";
    }
    
    double proposal_ratio = static_cast<double>(n_fi_sites) /
                            static_cast<double>(n_fu_sites);
    double P_acc = proposal_ratio * exp(-U);
    
    //---> accept or reject the move
    if(P_acc > thermal ) {
        //--- Accepted
        m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
        
        if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){   
            m_pState->GetVAHGlobalMeshProperties()->Add2Volume(new_Tvolume - old_Tvolume);
            m_pState->GetVAHGlobalMeshProperties()->Add2TotalArea(new_Tarea - old_Tarea);
            m_pState->GetVAHGlobalMeshProperties()->Add2GlobalCurvature(new_Tcurvature - old_Tcurvature);
        }
        return true;
    }
    else{ // reject the move

        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        for (auto it_l : Affected_links_old){
            it_l->Reverse_InteractionEnergy();
            it_l->Reverse_VFInteractionEnergy();
        }
        ReverseAScission(f_site);
        return false;
    }
    
    return false;
}
//========================================================================
//=====================  Fusion function =================================
//========================================================================
bool TopologyChangeByTriangularPrism::FusionByMove(fusion_site &fu_site, double thermal, int no_fus_sites){
    
    // There is a lot to be added to this function
       // if(m_AcceptedMoves!=0)
         //return false;

    double new_energy = 0;
    double old_energy = 0;
    
    triangle* t1 = fu_site.t1;
    triangle* t2 = fu_site.t2;
    
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
        if(!pV->SetCopy()){
            Nfunction::ConsolePrint_Error(" Error-> copying faild ");
        }        
        old_energy += pV->GetEnergy();
        old_energy += pV->GetBindingEnergy(); // vector field
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
            if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
                m_pState->GetVAHGlobalMeshProperties()->CalculateAPrismBasesContributionToGlobalVariables(t1, t2, old_Tvolume, old_Tarea, old_Tcurvature);
            }
        //====
           // we know perform fusion and our mesh is no longer the old one.
           // this will also updates trinagule and links normal and shape operator 

             fusion_outcome fusion_data;          
            if (!Fuse_MeshViaTwoTriangles(fu_site, fusion_data)) {
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
            
            
            if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
                m_pState->GetVAHGlobalMeshProperties()->CalculateAPrismFacesContributionToGlobalVariables(fusion_data.vertices, fusion_data.newLinks, new_Tvolume, new_Tarea, new_Tcurvature);
            }
               //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature); 
        
    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy + dE_volume + dE_t_area + dE_g_curv;

    
    tot_diff_energy -= m_Local_KG * 4 * PI;
    double U = m_Beta * tot_diff_energy - m_DBeta;
    
    std::vector<fission_site> fis_sites  = FindNecks();
    double no_fis_sites = fis_sites.size(); // this is M+1, because we get this is right after the move is done
    double proposal_ratio = static_cast<double>(no_fus_sites) /
                            static_cast<double>(no_fis_sites); // N/(M+1)
    
    double P_acc = proposal_ratio * exp(-U);
    
    //---> accept or reject the move
    if(P_acc > thermal ) {
            m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
            m_pState->GetVAHGlobalMeshProperties()->Add2Volume(new_Tvolume - old_Tvolume);
            m_pState->GetVAHGlobalMeshProperties()->Add2TotalArea(new_Tarea - old_Tarea);
            m_pState->GetVAHGlobalMeshProperties()->Add2GlobalCurvature(new_Tcurvature - old_Tcurvature);
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
std::vector<vertex*> V_vertices = {v1, v2, v3};

for (vertex* v : V_vertices) {
    std::vector<vertex*> neighbors = v->GetVNeighbourVertex();

    for (vertex* n : neighbors) {
        if (n == u1 || n == u2 || n == u3) {
            return false;
        }
    }
}
  
    return true;
        
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
//------------------------------------------------------------------------------
// Reverse_Fuse_MeshViaTwoTriangles
//------------------------------------------------------------------------------
// Reverts a previously successful Fuse_MeshViaTwoTriangles() operation.
//
// IMPORTANT:
//   - This function is intended exclusively for Monte Carlo rejection.
//   - It must be called  after a successful Fuse_MeshViaTwoTriangles() execution, before any mesh chanege 
//   - The provided fusion_outcome object must be the exact outcome returned
//     by that fusion operation.
//   - The function assumes that no other mesh modifications have occurred
//     since the fusion; otherwise the mesh state may become inconsistent.
//
// Restores:
//   - Original triangles (t1, t2)
//   - Original vertex and link states
//   - Ghost and active triangle/link pools
//   - Mirror-link containers
//------------------------------------------------------------------------------
    
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
bool TopologyChangeByTriangularPrism::Fuse_MeshViaTwoTriangles(fusion_site &fu_site, fusion_outcome &outcome){
//------------------------------------------------------------------------------
// Fuse_MeshViaTwoTriangles (June 2026, FreeDTS v2, Weria Pezeshkian)
//------------------------------------------------------------------------------
// Performs a local topological fusion of two triangles into a triangular-prism
// configuration.
//
// Workflow:
//   1. Verify sufficient ghost triangles and links are available.
//   2. Remove the two triangles from the active mesh.
//   3. Detach source triangles from their vertices.
//   4. Construct six new prism triangles using the supplied "topology map".
//   5. Reuse and update existing links; create new links from ghost pools.
//   6. Establish link connectivity and mirror-link relationships.
//   7. Update triangle normals, link normals, and link shape operators.
//   8. Populate the fusion outcome structure for potential rollback or analysis.
//
// Inputs:
//   fu_site   : Fusion site containing the two triangles and prism topology.
//
// Outputs:
//   outcome   : Stores all newly created and modified mesh entities.
//
// Returns:
//   true  - Fusion completed successfully.
//   false - Fusion failed due to insufficient ghost resources, invalid
//           topology, or mirror-link construction errors. false means error generally.
//
// Notes:
//   - Consumes 6 ghost triangles and 12 ghost links.
//   - Reuses existing edge links where possible.
//   - Assumes the supplied topology mapping is valid and consistent. No checking for it.
//------------------------------------------------------------------------------
         
if (m_pGhostT.size() < 6 || m_pGhostL.size() < 12) {
    std::cout << "---> Warning: Ghost counters are full. Fusion will not be performed. "
              << "This is a code limitation, not a physics limitation. "
              << "Please restart the simulation to restore correct physics behavior.\n";
    return false;
}
    
    triangle* t1 = fu_site.t1;
    triangle* t2 = fu_site.t2;
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

//--> remove the fusing triangles from the vertex list
    for(int i = 0; i<6; i++){
     if(i<3){ 
            Vver[i]->RemoveFromTraingleList(t1);
      } 
      else{
            Vver[i]->RemoveFromTraingleList(t2);
      }
    }
    std::vector<triple> VTriples = (fu_site.topology).VTriples;
    std::vector<triangle*> pnewTrinagles;       pnewTrinagles.reserve(6);
    std::vector<links*> pnewLinks;              pnewLinks.reserve(12);
    std::vector<links*> pHalfnewLinks;          pHalfnewLinks.reserve(6);
    std::vector<links*> poldLinks;              poldLinks.reserve(6);

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
                return false;
            }
            old_link->UpdateV3(tv3);  // very important 
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
            
            tv2->AddtoNeighbourVertex(tv3);
            tv3->AddtoNeighbourVertex(tv1);


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
            std::string message = "---> In " + TopologyChangeByTriangularPrism::GetDefaultReadName() + " command,  function Fuse_MeshViaTwoTrinagles() \n";
            Nfunction::ConsolePrint_Error(message);
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
    
    outcome.newTriangles = pnewTrinagles;
    outcome.newLinks = pnewLinks;
    outcome.oldLinks = poldLinks;
    outcome.pHalfnewLinks = pHalfnewLinks;
    outcome.vertices = Vver;
    outcome.t1 = t1;
    outcome.t2 = t2;

    return true;
}

// it creates a triangle and place it to the active trinagles list
triangle * TopologyChangeByTriangularPrism::CreateTriangleByTriple(
                    std::array<vertex*, 3> &v_ver, std::array<links*, 3> &v_links){

//---> create the triangle
    //--- select the last ghost triangle
    if(m_pGhostT.size() == 0){
        Nfunction::ConsolePrint_Error(" ---> error: (204847)\n");
    }
    triangle *gt1 = m_pGhostT[m_pGhostT.size()-1];
    //--- update its vertex
    gt1->UpdateVertex(v_ver[0], v_ver[1], v_ver[2]);

    //--- set the triangle from ghost to active
    m_pActiveT.push_back(gt1);
    m_pGhostT.pop_back();
    
//---> now make the changes in the edges
    //--- update their triangles
    for( auto it_l : v_links){
        it_l->UpdateTriangle(gt1);
    }
    for(int i = 0 ; i<3 ; i++){
        //---- update their next two edges
        v_links[i]->UpdateNeighborLink1(v_links[(i+1)%3]);
        v_links[i]->UpdateNeighborLink2(v_links[(i+2)%3]);
        //--- update their v3 (v2 and v1 remain the same)
        v_links[i]->UpdateV3(v_ver[(i+2)%3]);
    }    
//---> update the vertices triangle list
    for( auto it_v : v_ver){
        it_v->AddtoTraingleList(gt1);
    }
//---> update geometry; we cannot update the geometry of the the vertices because they have some edges that should be removed
    //--- update the normal and area of the triangle
    gt1->UpdateNormal_Area(&m_Box);   //trinagule normal and area should be obtained
    //--- update the three links notmal and shape Operator
    
    for( auto it_l : v_links){
        it_l->UpdateNormal();
        it_l->UpdateShapeOperator(&m_Box);

    }

    return gt1;
}
// this function cuts the neck made of p1 and p2, i.e., pair
std::vector <triangle *> TopologyChangeByTriangularPrism::DoAScission(fission_site &f_site){
    
    if(m_pGhostT.size()<2){
        std::cout<<" ---> [not an error] not enough reserved trinagles; you may restart the simulation "<<std::endl;
        exit(0);
    }  
    //--- create the triangles
    std::vector <triangle *> createdtriangles;
    createdtriangles.push_back(CreateTriangleByTriple(f_site.v_ver , f_site.V_links ));
    createdtriangles.push_back(CreateTriangleByTriple(f_site.u_ver , f_site.U_links));
    
    // removing the trinagles     
    for(auto it_l : f_site.C_Links){
        MakeATriangleGhost(it_l->GetTriangle());
    }
      // removing the edges        
    for(auto it_l : f_site.C_Links){
        MakeALinkGhost(it_l);
    }   
    for (auto it_v : f_site.v_ver){
        (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(it_v);
    }
    for (auto it_v : f_site.u_ver){
        (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(it_v);
    }
    
        return createdtriangles;
}
// this is the exact reverse action of DoAScission; different from DoAFussion
bool TopologyChangeByTriangularPrism::ReverseAScission(fission_site &f_site){

    for (auto it_v : f_site.v_ver){
        it_v->Reverse2PreviousCopy();
    }
    for (auto it_v : f_site.u_ver){
        it_v->Reverse2PreviousCopy();
    }
    for (auto it_l : f_site.V_links){
        it_l->Reverse2PreviousCopy();
        it_l->Reverse_VFInteractionEnergy();
    }    
    for (auto it_l : f_site.U_links){
        it_l->Reverse2PreviousCopy();
        it_l->Reverse_VFInteractionEnergy();
    }  
//-- add back the removed trinagles  
//-----> remove the added tringale   
    triangle * t1 = m_pActiveT.back();
    m_pActiveT.pop_back();
    triangle * t2 = m_pActiveT.back();
    m_pActiveT.pop_back();
//----> add the removed 6 trinagle   
    m_pActiveT.insert(m_pActiveT.end(),std::make_move_iterator(m_pGhostT.end() - 6),std::make_move_iterator(m_pGhostT.end()));
    m_pGhostT.resize(m_pGhostT.size() - 6);
//---- ghost both added trinagles
    m_pGhostT.push_back(t1);
    m_pGhostT.push_back(t2);
    
    
//--- return all the edges to life    
    m_pActiveL.insert(m_pActiveL.end(),std::make_move_iterator(m_pGhostL.end() - 12),std::make_move_iterator(m_pGhostL.end()));
    m_pGhostL.resize(m_pGhostL.size() - 12);
    
    for (auto it_l : f_site.C_Links){
        auto it_m = it_l->GetMirrorLink();
        m_pLeftL.push_back(it_l);
        m_pRightL.push_back(it_m);
    }         
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
auto canonicalize_vertices = [](vertex*& a, vertex*& b, vertex*& c) {
    if (a > b) std::swap(a, b);
    if (a > c) std::swap(a, c);
    if (b > c) std::swap(b, c);
};
//== this function get a mesh and search through it and finds possible fission sites
std::vector<fission_site> TopologyChangeByTriangularPrism::FindNecks() {
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
    std::vector<fission_site> pair_list;

    std::unordered_set<PotTriangleKey, PotTriangleKeyHash> unique_tris;

    int id = 0;

    for (auto *pv : m_pSurfV) {
        

        auto nv = pv->GetVNeighbourVertex();

        for (size_t i = 0; i < nv.size(); i++) {
            vertex* pv2 = nv[i];
            for (size_t j = i + 1; j < nv.size(); j++) {
                vertex* pv3 = nv[j];
                vertex* pv1 = pv;
                if (!pv1->IsAnyTriangle(pv2, pv3)) {  // vertices are nighbour of pv1 but not forming a triangle already
                    if(!pv2->IsThereAConnectingLink(pv3)) // pv2-pv3 are connected
                            continue;
                    

                   // ---- canonical ordering
                    canonicalize_vertices(pv1, pv2, pv3);  // based on pointer value not location
                        if(!(pv1->IsThereAConnectingLink(pv2)) ){
        Nfunction::ConsolePrint_Error("--> error what the hell \n");
        exit(0);
    }
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

    for (int i = 0; i < (int)trinagle_frame.size(); i++) {
        for (int j = i + 1; j < (int)trinagle_frame.size(); j++) {

            fission_site temp;

            if (Is_A_Neck(trinagle_frame[i], trinagle_frame[j], temp)) {
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
bool TopologyChangeByTriangularPrism::Is_A_Neck(pot_triangle potT1, pot_triangle potT2, fission_site & prism_site) {
    /**
     * @brief Checks if two potential triangles (pot_triangle) are well-connected to form a neck for fission.
*/

std::array<vertex*, 3> v_ver = {potT1.pv1, potT1.pv2, potT1.pv3};
std::array<vertex*, 3> u_ver = {potT2.pv1, potT2.pv2, potT2.pv3};

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
                connections[i][j] = 1;
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
// the system is not a prism
if(bride_links.size() != 6){
   return false;
}
// now fix the orientation and get the links
    links *lv1_2;
    
    if(!v_ver[0]->IsThereAConnectingLink(v_ver[1], lv1_2)){
        Nfunction::ConsolePrint_Error("--> error (v_ver0 is not connected to v_ver1): this should not happen \n");
    }

    if(lv1_2->GetV3() != u_ver[0] && lv1_2->GetV3() != u_ver[1] && lv1_2->GetV3() != u_ver[2]){
         v_ver[0] = potT1.pv2;   
         v_ver[1] = potT1.pv1; 
         if(lv1_2->GetMirrorFlag()){
            lv1_2 = lv1_2->GetMirrorLink();
        }
        else{
            Nfunction::ConsolePrint_Error("--> error (47473): this should not be tested \n");
        }
    }
    
    links *lu1_2;
    if(!u_ver[0]->IsThereAConnectingLink(u_ver[1], lu1_2)){
        Nfunction::ConsolePrint_Error("--> error (47744): this should not happen \n");
    }

    if(lu1_2->GetV3() != v_ver[0] && lu1_2->GetV3() != v_ver[1] && lu1_2->GetV3() != v_ver[2]){
         u_ver[0] = potT2.pv2;   
         u_ver[1] = potT2.pv1; 
         if(lu1_2->GetMirrorFlag()){
            lu1_2 = lu1_2->GetMirrorLink();
        }
        else{
            Nfunction::ConsolePrint_Error("--> error (47473): this should not be tested \n");
        }
    }
 //== now orintation is fixed the below links does not need orintation fix because the trinagle is already fixed.   
      links *lv2_3 = v_ver[1]->GetConnectingLink(v_ver[2]);
      links *lv3_1 = v_ver[2]->GetConnectingLink(v_ver[0]);
      links *lu2_3 = u_ver[1]->GetConnectingLink(u_ver[2]);
      links *lu3_1 = u_ver[2]->GetConnectingLink(u_ver[0]);
      
        
    std::array<links*, 3> v_links = {lv1_2, lv2_3, lv3_1};
    std::array<links*, 3> u_links = {lu1_2, lu2_3, lu3_1};
    
    
    prism_site.v_ver = v_ver;
    prism_site.u_ver = u_ver;
    prism_site.C_Links = bride_links;
    prism_site.V_links = v_links;
    prism_site.U_links = u_links;
    
    if(!CheckFaceAngleForFissionSite(prism_site)){
        return false;
    }
    return true;


}
bool TopologyChangeByTriangularPrism::CheckFaceAngleForFissionSite(fission_site &f_site){
    
    triangle t1(-1, f_site.v_ver[0], f_site.v_ver[1], f_site.v_ver[2]);
    triangle t2(-2, f_site.u_ver[0], f_site.u_ver[1], f_site.u_ver[2]);
    Vec3D n1 = t1.CalculateNormal(m_Box);
    Vec3D n2 = t2.CalculateNormal(m_Box);

    for(auto it_l : f_site.V_links){
        Vec3D normal_v = it_l->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n1,normal_v) < m_MinAngle ){
            return false;
        }
    }
    
    for(auto it_l : f_site.U_links){
        Vec3D normal_v = it_l->GetMirrorLink()->GetTriangle()->GetNormalVector();
        if( n1.dot(n2,normal_v) < m_MinAngle ){
            return false;
        }
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
std::vector<links*> TopologyChangeByTriangularPrism::GetEdgesWithInteractionChange(fission_site &f_site){
/**
 * @brief Collects links whose interaction energies may change during a scission move.
 *
 * This function identifies the subset of mesh links whose inclusion or vector-field
 * interaction energies must be recomputed after a topology modification at the
 * specified fission site.
 *
 * The search is restricted to links connected to vertices participating in the
 * scission (`f_site.v_ver` and `f_site.u_ver`). A link is considered potentially
 * affected if:
 * - The associated vertex hosts an inclusion, or
 * - The system contains one or more vector-field layers.
 *
 * All neighboring links of such vertices are gathered and then filtered to
 * produce a unique set of links. Duplicate entries and mirror-link pairs are
 * removed to ensure that each physical interaction is represented only once.
 *
 * The returned links are typically used before and after a trial topology
 * change to:
 * - Store current interaction energies,
 * - Recompute inclusion interaction energies,
 * - Recompute vector-field interaction energies, and
 * - Restore energies if the move is rejected.
 *
 * @param f_site
 *        Description of the candidate scission region, including the vertices
 *        whose local connectivity and interactions may be modified.
 *
 * @return A unique collection of links whose interaction energy contributions
 *         may be affected by the proposed topology change.
 *
 * @note In development mode, a diagnostic message is printed indicating that
 *       additional pre-filtering optimizations may be possible.
 */
#if DEVELOPMENT_MODE == Enabled
    std::cout<<" DEVELOPMENT_MODE ID 9595473: This function can be  made  better\n ";
    // for example, precheck to select only links that both has includioon ...
#endif
    
    std::vector<links*> edge_with_interaction_change;
    std::vector<links *> all_temlinks;

    bool system_has_vf = false;
    if(m_No_Vectorfield != 0){
        system_has_vf = true;
    }
    
    for (auto v_it : f_site.v_ver){
        if(v_it->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v_it->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
        }
    }  
    for (auto v_it : f_site.u_ver){
        if(v_it->VertexOwnInclusion() || system_has_vf ) {  // due to vector fields
        std::vector <links *> n_edges = v_it->GetVLinkList();
        all_temlinks.insert(all_temlinks.end(), n_edges.begin(), n_edges.end());
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
    state += " "+ Nfunction::D2S(m_Period) + " "+ m_PrismMapTopologyFile;
    return state;
}
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
