

#include <stdio.h>
#include "PositionRescaleFrameTensionCoupling.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"


/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
PositionRescaleFrameTensionCoupling::PositionRescaleFrameTensionCoupling(int period, double sigma, std::string direction, State *pState)  :
        m_pState(pState),
        m_pActiveV(pState->GetMesh()->GetActiveV()),
        m_pActiveT(pState->GetMesh()->GetActiveT()),
        m_pRightL(pState->GetMesh()->GetRightL()),
        m_pEdgeL(pState->GetMesh()->GetEdgeL()),
        m_Beta(pState->GetSimulation()->GetBeta()),
        m_DBeta(pState->GetSimulation()->GetDBeta()),
        m_MinLength2(pState->GetSimulation()->GetMinL2()),
        m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
        m_MinAngle(pState->GetSimulation()->GetMinAngle()),
        m_Type (direction){
            
            m_SigmaP = sigma;
            m_Period = period;
            //-- convert the type string into the direction vector
            SetDirection(m_Type);

}
PositionRescaleFrameTensionCoupling::~PositionRescaleFrameTensionCoupling() {
    
}
void PositionRescaleFrameTensionCoupling::Initialize() {
    
    std::cout<<"---> the algorithm for box change involves applying: "<< GetDefaultReadName()<<" \n";
    m_pBox = (m_pState->GetMesh())->GetBox();
}
bool PositionRescaleFrameTensionCoupling::ChangeBoxSize(int step){
    /**
     * @brief a call function to change the simulation box size at a given step.
     *
     * This function changes the size of the simulation box based on the current step.
     * It first checks if the current step matches the defined period. If the voxel size
     * is below a threshold, it updates the voxel size and re-voxelizes. It then computes
     * the size change for the box using, and attempts to change the
     * box size based on these calculations.
     *
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
//---> if does not match the preiod, return false
    if(step%m_Period != 0)
        return false;
    //--- first check if the voxel size is fine
    double voxel_lx = m_pState->GetVoxelization()->GetXSideVoxel();
    double voxel_ly = m_pState->GetVoxelization()->GetYSideVoxel();
    double voxel_lz = m_pState->GetVoxelization()->GetZSideVoxel();
    if(voxel_lx < 1.05 || voxel_ly < 1.05 || voxel_lz < 1.05) {
        m_pState->GetVoxelization()->UpdateVoxelSize(1.2, 1.2, 1.2);
        m_pState->GetVoxelization()->Voxelize(m_pActiveV);
    }
    
//---> find the size of box change; isotropic method (here all the other methods can be performed)
    double dx = 1 - 2 * (m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
    dx *= m_DR;
    double dy = dx*((*m_pBox)(1))/(*m_pBox)(0);
    double dz = dx*((*m_pBox)(2))/(*m_pBox)(0);
    
    // lx, ly, lz how much the box should be scaled in each direction
     double lx = 1 + m_Direction(0) * dx / (*m_pBox)(0);
     double ly = 1 + m_Direction(1) * dy / (*m_pBox)(1);
     double lz = 1 + m_Direction(2) * dz / (*m_pBox)(2);
    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    
    // Attempt to change the box size
    if (AnAtemptToChangeBox(lx, ly, lz, thermal)) {
        m_AcceptedMoves++;
        m_NumberOfAttemptedMoves++;
    }
    else{
        m_NumberOfAttemptedMoves++;
    }
    
    return true;
}
bool PositionRescaleFrameTensionCoupling::AnAtemptToChangeBox(double lx,double ly, double lz, double temp){

//---> check if we do the move, the distance will be normal
    if(!VertexMoveIsFine(lx, ly, lz)){
        return false;
    }

    //---> get the energies
    double old_energy = m_pState->GetEnergyCalculator()->GetEnergy();
    double new_energy = 0;
    
    double new_systemsize = 1;   // area for 2d, volume for 3d, line for 1d
    double old_systemsize = 1;
    for (int i=0;i<3;i++) {
        if(m_Direction(i) != 0){
            old_systemsize *= (*m_pBox)(i);
        }
    }
    
    // Obtain and sum the initial global variables that might change
    double old_Tvolume = 0.0, old_Tarea = 0.0, old_Tcurvature = 0.0;
    double new_Tvolume = 0.0, new_Tarea = 0.0, new_Tcurvature = 0.0;
//--->
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){

        old_Tvolume =  m_pState->GetVAHGlobalMeshProperties()->GetTotalVolume();
        old_Tarea =  m_pState->GetVAHGlobalMeshProperties()->GetTotalArea();
        old_Tcurvature =  m_pState->GetVAHGlobalMeshProperties()->GetTotalMeanCurvature();
    }
    for (std::vector<triangle *>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it){
        (*it)->ConstantMesh_Copy();
    }
    for (std::vector<links *>::iterator it = m_pRightL.begin() ; it != m_pRightL.end(); ++it){
        (*it)->ConstantMesh_Copy();
        (*it)->Copy_VFInteractionEnergy();
    }
    for (std::vector<links *>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it){
        (*it)->ConstantMesh_Copy();
        (*it)->Copy_VFInteractionEnergy();
    }
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        (*it)->ConstantMesh_Copy();
        (*it)->Copy_VFsBindingEnergy();
    }
    //---> for now, only active nematic force: ForceonVerticesfromInclusions
    // if we do this here, the force is from inital configuration
    // while if we do it after the move, it will be from final configurations
    double dE_force_from_inc = 0;
    double dE_force_from_vectorfields = 0;

    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        double x = (*it)->GetVXPos();
        double y = (*it)->GetVYPos();
        double z = (*it)->GetVYPos();
        double dx = x * (lx - 1);
        double dy = y * (ly - 1);
        double dz = z * (lz - 1);
        
        Vec3D DX(dx,dy, dz);
        dE_force_from_inc  += m_pState->GetForceonVerticesfromInclusions()->Energy_of_Force( *it, DX);
        dE_force_from_vectorfields += m_pState->GetForceonVerticesfromVectorFields()->Energy_of_Force( *it, DX);
    }
    
    (*m_pBox)(0) *= lx;
    (*m_pBox)(1) *= ly;
    (*m_pBox)(2) *= lz;
    for (std::vector<vertex*>::iterator it =  m_pActiveV.begin(); it != m_pActiveV.end(); ++it) {
        (*it)->ScalePos(lx,ly,lz);
    }
    //-- Update geometry
    for (std::vector<triangle *>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it){
        (*it)->UpdateNormal_Area(m_pBox);
    }

    if(!CheckFaces()){
        for (std::vector<vertex*>::iterator it =  m_pActiveV.begin(); it != m_pActiveV.end(); ++it) {
            (*it)->ScalePos(1.0/lx,1.0/ly,1.0/lz);
        }
        for (std::vector<triangle *>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        return false;
    }
    //--- lets now do the move
  //  m_pState->GetCurvatureCalculator()->Initialize();
    for (std::vector<links *>::iterator it = m_pRightL.begin() ; it != m_pRightL.end(); ++it){
        (*it)->UpdateShapeOperator(m_pBox);
    }
    for (std::vector<links *>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it){
        (*it)->UpdateEdgeVector(m_pBox);  //
    }
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(*it);
    }
   // instead of the we do below. so for paralization helps
 //new_energy = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();

//--- calculate new energies
    int number_of_vectorfields = m_pState->GetMesh()->GetNoVFPerVertex();
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it) {
        new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(*it);
        new_energy += (m_pState->GetEnergyCalculator())->CalculateVectorFieldMembraneBindingEnergy(*it);
    }
    for (std::vector<links *>::iterator it = m_pRightL.begin() ; it != m_pRightL.end(); ++it) {
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        for (int i = 0; i<number_of_vectorfields; i++){
            new_energy += (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(i, *it);
        }
    }
    for (std::vector<links *>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it) {
        new_energy += (m_pState->GetEnergyCalculator())->TwoInclusionsInteractionEnergy(*it);
        for (int i = 0; i<number_of_vectorfields; i++){
            new_energy += (m_pState->GetEnergyCalculator())->TwoVectorFieldInteractionEnergy(i, *it);
        }
    }
    
    //-- more to energies
    //---> get energy for ApplyConstraintBetweenGroups
    double dE_Cgroup = m_pState->GetApplyConstraintBetweenGroups()->CalculateEnergyChange(lx, ly, lz);


//---> new global variables
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateBoxRescalingContributionToGlobalVariables(lx, ly, lz, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    
    //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
    
    //--> only elatsic energy
    double diff_energy = new_energy - old_energy;
    //--> sum of all the energies
    double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_from_inc + dE_volume + dE_t_area + dE_g_curv ;
    double NV = m_pActiveV.size();
    
    for (int i=0;i<3;i++) {
        if(m_Direction(i) != 0){
            new_systemsize *= (*m_pBox)(i);
        }
    }
    
    tot_diff_energy -= m_SigmaP * (new_systemsize - old_systemsize);
    


    //---> accept or reject the move

    if ( pow(lx*ly*lz , NV) * exp(-m_Beta * tot_diff_energy + m_DBeta) > temp ) {
        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        
        //---> ApplyConstraintBetweenGroups
         m_pState->GetApplyConstraintBetweenGroups()->AcceptMove();
        
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
        //---> reverse the triangles
        (*m_pBox)(0) *= 1/lx;
        (*m_pBox)(1) *= 1/ly;
        (*m_pBox)(2) *= 1/lz;
        for (std::vector<triangle *>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
        }
        for (std::vector<links *>::iterator it = m_pRightL.begin() ; it != m_pRightL.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
            (*it)->Reverse_VFInteractionEnergy();
        }
        for (std::vector<links *>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
            (*it)->Reverse_VFInteractionEnergy();
        }
        for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
            (*it)->ReverseConstantMesh_Copy();
            (*it)->Reverse_VFsBindingEnergy();
        }
        return false;
     }

    return true;
}
bool PositionRescaleFrameTensionCoupling::VertexMoveIsFine(double lx,double ly, double lz){
    /**
     * @brief Checks if distances are valid given the new box dimensions.
     *
     * and ensuring the distances between vertices are within acceptable limits.
     *
     * @param lx scaling factor in the x-direction.
     * @param ly scaling factor in the y-direction.
     * @param lz scaling factor in the z-direction.
     * @return true if it is good
     */
    if(!CheckLinkLength(lx, ly, lz)){
        return false;
    }
    
    //--- if the Stretching is positive, so the distances should be fine
    if(lx>=1 && ly>=1 && lz>=1)
        return true;
    
    //--- checking the distance between each pair
    Voxel<vertex>  ****voxels = m_pState->GetVoxelization()->GetAllVoxel();
    int Nx = m_pState->GetVoxelization()->GetXVoxelNumber();
    int Ny = m_pState->GetVoxelization()->GetYVoxelNumber();
    int Nz = m_pState->GetVoxelization()->GetZVoxelNumber();

    for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
    for (int k = 0; k < Nz; ++k) {
        
        std::vector<vertex*> voxel_ver = (voxels[i][j][k])->GetContentObjects();
            
        //--- check distances within the same voxel
        for (std::vector<vertex*>::iterator it1 = voxel_ver.begin(); it1 != voxel_ver.end(); ++it1) {
            for (std::vector<vertex*>::iterator it2 = it1 + 1; it2 != voxel_ver.end(); ++it2) {
                    double l2 = StretchedDistanceSquardBetweenTwoVertices(*it1, *it2, lx, ly, lz);
                    if (l2 < m_MinLength2) {
                        return false;
                    }
                }
        }
        // check distances of vertex from neighbouring voxels
            for (int n = 0; n < 2; ++n)
            for (int m = 0; m < 2; ++m)
            for (int s = 0; s < 2; ++s)
              if(n != 0 || m != 0 || s!=0) {
                    std::vector<vertex*> voxel_ver2 = (voxels[i][j][k])->GetANeighbourCell(n,m,s)->GetContentObjects();

                    for (std::vector<vertex *>::iterator it1 = voxel_ver.begin() ; it1 != voxel_ver.end(); ++it1) {
                        for (std::vector<vertex *>::iterator it2 = voxel_ver2.begin() ; it2 != voxel_ver2.end(); ++it2) {
                            if(it1 != it2){
                                double l2 = StretchedDistanceSquardBetweenTwoVertices(*it1, *it2, lx, ly, lz);
                                if (l2 < m_MinLength2) {
                                    return false;
                                }
                        }//  if(it1 != it2){
                    }
                }
            }// for (int s = 0; s < 2; ++s) {
    }
    }
    }
    
    return true;
}
bool PositionRescaleFrameTensionCoupling::CheckLinkLength(double lx,double ly, double lz){
    /**
     * Checks if the length of links in the mesh is within acceptable bounds after rescaling.
     *
     * This function iterates through all links in the edge and right link lists of the mesh,
     * calculating their squared distances after rescaling. If any link's length squared is
     * outside the acceptable range, the function returns false.
     *
     * @param lx Scaling factor for the x-axis.
     * @param ly Scaling factor for the y-axis.
     * @param lz Scaling factor for the z-axis.
     * @return True if all links are within the acceptable length range, false otherwise.
     */
    
    for (std::vector<links *>::iterator it = m_pEdgeL.begin(); it != m_pEdgeL.end(); ++it) {
        double l2 = StretchedDistanceSquardBetweenTwoVertices((*it)->GetV1(), (*it)->GetV2(), lx, ly, lz);
        if (l2 < m_MinLength2 || l2 > m_MaxLength2) {
            return false;
        }
    }

    for (std::vector<links *>::iterator it = m_pRightL.begin(); it != m_pRightL.end(); ++it) {
        double l2 = StretchedDistanceSquardBetweenTwoVertices((*it)->GetV1(), (*it)->GetV2(), lx, ly, lz);
        if (l2 < m_MinLength2 || l2 > m_MaxLength2) {
            return false;
        }
    }

    return true;
}
double PositionRescaleFrameTensionCoupling::StretchedDistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2, double lx, double ly, double lz){
    /**
     * Calculates the squared stretched distance between two vertices.
     *
     * This function computes the squared distance between two vertices (v1 and v2)
     * after applying a stretching factor along each axis (lx, ly, lz). The stretching
     * factors allow for scaling the distance calculation to account for non-uniform
     * scaling of the simulation box dimensions. The function also ensures that the
     * distance respects periodic boundary conditions, wrapping around the box if
     * necessary.
     *
     * @param v1 Pointer to the first vertex.
     * @param v2 Pointer to the second vertex.
     * @param lx Scaling factor for the x-axis.
     * @param ly Scaling factor for the y-axis.
     * @param lz Scaling factor for the z-axis.
     * @return The squared distance between the two vertices after scaling and
     *         considering periodic boundary conditions.
     */
    // Calculate the initial distance components between the two vertices
    double dx = v2->GetVXPos() - v1->GetVXPos();
    double dy = v2->GetVYPos() - v1->GetVYPos();
    double dz = v2->GetVZPos() - v1->GetVZPos();
    
    // Apply the stretching factors to each component
    dx *= lx;
    dy *= ly;
    dz *= lz;

    // Calculate the dimensions of the scaled simulation box
    double Lx = lx * (*m_pBox)(0);
    double Ly = ly * (*m_pBox)(1);
    double Lz = lz * (*m_pBox)(2);

    // Apply periodic boundary conditions for the x-axis
    if(fabs(dx) > Lx / 2.0) {
        if(dx < 0)
            dx = Lx + dx;
        else
            dx = dx - Lx;
    }
    // Apply periodic boundary conditions for the y-axis
    if(fabs(dy) > Ly / 2.0) {
        if(dy < 0)
            dy = Ly + dy;
        else
            dy = dy - Ly;
    }
    // Apply periodic boundary conditions for the z-axis
    if(fabs(dz) > Lz / 2.0) {
        if(dz < 0)
            dz = Lz + dz;
        else
            dz = dz - Lz;
    }
    // Return the squared distance after scaling and applying periodic boundary conditions
    return dx * dx + dy * dy + dz * dz;
}
bool PositionRescaleFrameTensionCoupling::CheckFaceAngleOfOneLink(links* p_edge) {
    
    // Function to check if the angle between the normal vectors of the faces sharing the given edge
    // is above a minimum threshold defined by m_MinAngle.
    
    if(!p_edge->GetMirrorFlag()){
        return true;
    }
    // Retrieve the normal vector of the triangle associated with the edge.
    Vec3D n1 = p_edge->GetTriangle()->GetNormalVector();

    // Retrieve the normal vector of the triangle associated with the mirror edge.
   
    Vec3D n2 = p_edge->GetMirrorLink()->GetTriangle()->GetNormalVector();

    // Check if the dot product of the two normal vectors is less than the minimum allowed angle.
    // If it is, return false indicating that the face angle condition is not satisfied.
    if (n1.dot(n1, n2) < m_MinAngle) {
        return false;
    }

    // If the angle condition is satisfied, return true.
    return true;
}

bool PositionRescaleFrameTensionCoupling::CheckFaces() {
    // Function to check if the angle between the faces of all links in the m_pRightL list
    // satisfies the minimum angle condition.
    
    // Iterate through all links in the m_pRightL vector.
    for (std::vector<links*>::iterator it = m_pRightL.begin(); it != m_pRightL.end(); ++it) {
        // For each link, check if the face angle condition is satisfied using CheckFaceAngleOfOneLink function.
        // If any link does not satisfy the condition, return false.
        if (!CheckFaceAngleOfOneLink(*it)) {
            return false;
        }
    }

    // If all links satisfy the face angle condition, return true.
    return true;
}
void PositionRescaleFrameTensionCoupling::SetDirection(std::string direction){
    /**
     * @brief Sets the direction for box size rescaling.
     *
     * This function sets the direction vector for rescaling the box size based on the
     * provided string. It supports various combinations of X, Y, and Z directions.
     *
     * @param direction A string specifying the direction for rescaling ("XYZ", "XY", "XZ", "YZ", "X", "Y", "Z").
     */
    if(direction == "XYZ"){
        m_Direction(0) = 1;
        m_Direction(1) = 1;
        m_Direction(2) = 1;

    }
    else if(direction == "XY"){
        m_Direction(0) = 1;
        m_Direction(1) = 1;
        m_Direction(2) = 0;

    }
    else if(direction == "XZ"){
        m_Direction(0) = 1;
        m_Direction(1) = 0;
        m_Direction(2) = 1;

    }
    else if(direction == "YZ"){
        m_Direction(0) = 0;
        m_Direction(1) = 1;
        m_Direction(2) = 1;

    }
    else if(direction == "X"){
        m_Direction(0) = 1;
        m_Direction(1) = 0;
        m_Direction(2) = 0;
    }
    else if(direction == "Y"){
        m_Direction(0) = 0;
        m_Direction(1) = 1;
        m_Direction(2) = 0;
    }
    else if(direction == "Z"){
        m_Direction(0) = 0;
        m_Direction(1) = 0;
        m_Direction(2) = 1;
    }
    else{
        // Print an error message if the direction is unknown
        std::cout << "---> Error: direction of the box change is unknown\n";    }
    
    return;
}
std::string PositionRescaleFrameTensionCoupling::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}




