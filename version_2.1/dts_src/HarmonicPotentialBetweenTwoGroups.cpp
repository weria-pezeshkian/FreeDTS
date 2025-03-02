


#include "HarmonicPotentialBetweenTwoGroups.h"
#include "Nfunction.h"
#include "State.h"
HarmonicPotentialBetweenTwoGroups::HarmonicPotentialBetweenTwoGroups(State* pState, double K, double l0, std::string group1,std::string group2,double nx,double ny,double nz) {
    m_K = K/2;
    m_Group1Name = group1;
    m_Group2Name = group2;
    m_Direction(0) = nx;
    m_Direction(1) = ny;
    m_Direction(2) = nz;
    m_L0 = l0;
    m_Energy = 0;
    m_G1Size = 1;
    m_G2Size = 1;
    m_DE = 0;
    m_pState = pState;
}
HarmonicPotentialBetweenTwoGroups::~HarmonicPotentialBetweenTwoGroups(){
    
}
bool HarmonicPotentialBetweenTwoGroups::Initialize() {
    /*
     * @brief Initializes the HarmonicPotentialBetweenTwoGroups object by setting up the groups and calculating initial parameters.
     *
     * This function initializes the object by retrieving the groups from the provided group names, calculating their centers of geomtry,
     * checking the distance between the groups, normalizing the direction vector, and calculating the initial energy.
     *
     * @return True if initialization is successful, false otherwise.
     */
    
    
    std::map<std::string, std::vector<vertex*> >& L_Groups = m_pState->GetMesh()->GetGroups();
    m_pBox = m_pState->GetMesh()->GetBox();
    // Check if both groups exist in the provided group map
    if (L_Groups.find(m_Group1Name) != L_Groups.end() && L_Groups.find(m_Group2Name) != L_Groups.end()) {
        // Retrieve the groups from the map
        m_pGroup1 = L_Groups.at(m_Group1Name);
        m_pGroup2 = L_Groups.at(m_Group2Name);
    } else {
        std::cout << "---error--> groups for " << GetDefaultReadName() << " do not exist.\n";
        return false;
    }

    // Calculate the center of geometry (COG) for both groups
    m_Group1COG = COMVertexGroup(m_pGroup1);
    m_Group2COG = COMVertexGroup(m_pGroup2);

    // Set the size of each group
    m_G1Size = static_cast<double>(m_pGroup1.size());
    m_G2Size = static_cast<double>(m_pGroup2.size());

    // Check if the distance between the groups is valid
    if (!DistanceCheck()) {
        return false;
    }

    // Normalize the direction vector
    m_Direction.normalize();

    // Calculate the initial distance and energy
    double dist = Vec3D::dot(m_Direction, m_Group2COG - m_Group1COG);
    m_Energy = m_K * (dist - m_L0) * (dist - m_L0);

    return true;
}
std::string HarmonicPotentialBetweenTwoGroups::Output_DataString(){

    double dist = Vec3D::dot(m_Direction, m_Group2COG - m_Group1COG);
    double force = -2*m_K * (dist - m_L0);
    std::string data_output = Nfunction::D2S(m_Energy) + " " + Nfunction::D2S(force)+ "  " + Nfunction::D2S(dist);
    
    return data_output;
}
bool HarmonicPotentialBetweenTwoGroups::DistanceCheck() {
    /*
     * @brief Checks if the distance between the centers of geomtry (COG) of two groups exceeds half the box size in any specified direction.
     *
     * This function computes the distance between the centers of geomtry of two groups and checks if it exceeds half of the box size in the specified directions (X, Y, Z).
     * If the distance exceeds half the box size in any of these directions, the simulation needs to be stopped.
     *
     * @return True if the distance is within the allowed range, false otherwise.
     */
    // Calculate the distance between the centers of gravity of the two groups
    Vec3D Dist = m_Group1COG - m_Group2COG;

    // Check if the distance in the X direction exceeds half of the box size
    if (m_Direction(0) != 0 && fabs(Dist(0)) > (*m_pBox)(0) / 2.0) {
        std::cout << "---> The pulling distance is larger than half of the box size in the X direction, simulation needs to be stopped." << std::endl;
        
        *(m_pState->GetTimeSeriesLog()) << "---> The pulling distance is larger than half of the box size in the X direction, simulation needs to be stopped. \n";

        return false;
    }
    // Check if the distance in the Y direction exceeds half of the box size
    else if (m_Direction(1) != 0 && fabs(Dist(1)) > (*m_pBox)(1) / 2.0) {
        std::cout << "---> The pulling distance is larger than half of the box size in the Y direction, simulation needs to be stopped." << std::endl;
        *(m_pState->GetTimeSeriesLog()) << "---> The pulling distance is larger than half of the box size in the Y direction, simulation needs to be stopped. \n";
        return false;
    }
    // Check if the distance in the Z direction exceeds half of the box size
    else if (m_Direction(2) != 0 && fabs(Dist(2)) > (*m_pBox)(2) / 2.0) {
        std::cout << "---> The pulling distance is larger than half of the box size in the Z direction, simulation needs to be stopped." << std::endl;
        *(m_pState->GetTimeSeriesLog()) << "---> The pulling distance is larger than half of the box size in the Z direction, simulation needs to be stopped. \n";
        return false;
    }

    // If none of the distances exceed the allowed limit, return true
    return true;
}
double HarmonicPotentialBetweenTwoGroups::CalculateEnergyChange(vertex* p_vertex, Vec3D Dx) {

    m_T_Group1COG = m_Group1COG;
    m_T_Group2COG = m_Group2COG;
    if(p_vertex->GetGroupName() == m_Group1Name) {
        
        m_T_Group1COG = m_T_Group1COG + Dx*(1.0/m_G1Size);
    }
    else if(p_vertex->GetGroupName() == m_Group2Name){
        
        m_T_Group2COG = m_T_Group2COG + Dx*(1.0/m_G2Size);
    }
    else { // meaning the vertex move does not affect the distance between the two group
        return 0;
    }
    if(!DistanceCheck()){
        exit(EXIT_FAILURE);
    }
    
    double dist = Vec3D::dot(m_Direction,m_T_Group2COG - m_T_Group1COG);
    m_DE = m_K*(dist-m_L0)*(dist-m_L0) - m_Energy;
    
    return m_DE;
}
double HarmonicPotentialBetweenTwoGroups::CalculateEnergyChange(double lx, double ly, double lz) {

    m_T_Group1COG(0) = m_Group1COG(0)*(lx);
    m_T_Group1COG(1) = m_Group1COG(1)*(ly);
    m_T_Group1COG(2) = m_Group1COG(2)*(lz);

    m_T_Group2COG(0) = m_Group2COG(0)*(lx);
    m_T_Group2COG(1) = m_Group2COG(1)*(ly);
    m_T_Group2COG(2) = m_Group2COG(2)*(lz);
    
    if(!DistanceCheck()){
        exit(EXIT_FAILURE);
    }
    
    double dist = Vec3D::dot(m_Direction,m_T_Group2COG - m_T_Group1COG);
    m_DE = m_K*(dist-m_L0)*(dist-m_L0) - m_Energy;
    
    return m_DE;
}
void HarmonicPotentialBetweenTwoGroups::AcceptMove()
{
    m_Group1COG = m_T_Group1COG;
    m_Group2COG = m_T_Group2COG;
    m_Energy += m_DE;
    return;
}
Vec3D HarmonicPotentialBetweenTwoGroups::COMVertexGroup(std::vector<vertex *> ver)
{
    Vec3D com;
    double x=0;
    double y=0;
    double z=0;
    for (std::vector<vertex *>::iterator it = ver.begin() ; it != ver.end(); ++it)
    {
        x+=(*it)->GetVXPos();
        y+=(*it)->GetVYPos();
        z+=(*it)->GetVZPos();
    }
    com(0)=x/ver.size();
    com(1)=y/ver.size();
    com(2)=z/ver.size();

    return com;
}



std::string HarmonicPotentialBetweenTwoGroups::CurrentState(){
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName() + " "+ Nfunction::D2S(2*m_K)+" "+Nfunction::D2S(m_L0);
        state = state +" "+ m_Group1Name +" "+ m_Group2Name +" "+Nfunction::D2S(m_Direction(0))+" "+Nfunction::D2S(m_Direction(1));
        state = state +" "+Nfunction::D2S(m_Direction(2));
        return state;
}

