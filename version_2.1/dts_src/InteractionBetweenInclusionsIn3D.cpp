#include <unordered_set>
#include "InteractionBetweenInclusionsIn3D.h"
#include "State.h"

// Constructor
InteractionBetweenInclusionsIn3D::InteractionBetweenInclusionsIn3D(State* pstate, std::string input_data)  : m_pState(pstate), m_Input_Data(input_data) {

}
InteractionBetweenInclusionsIn3D::~InteractionBetweenInclusionsIn3D() {
    // Destructor
}
void InteractionBetweenInclusionsIn3D::Initialize() {
    std::vector<std::string> data = Nfunction::Split(m_Input_Data);
    
    // Prevent out-of-bounds access
    if (data.size() < 1) {
        std::cerr << "---> error: insufficient input data for "
                  << this->GetDefaultReadName() << ".\n";
        exit(-1);
    }

    // Use correct parsing function (double instead of int)
    m_EP = Nfunction::String_to_Double(data[0]);


    // Improved error messages and checks
    if (m_EP < 0) {
        std::cout << "---> error: assigned interaction strength in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }

    double vx = m_pState->GetVoxelization()->GetXSideVoxel();
    double vy = m_pState->GetVoxelization()->GetYSideVoxel();
    double vz = m_pState->GetVoxelization()->GetZSideVoxel();
    
    //std::cout<<vx<<" "<<

    return;
}
double InteractionBetweenInclusionsIn3D::GetVertexNonBondedEnergy(vertex* pvertex) const {
    

    // Check if the vertex owns an inclusion
    if (!pvertex->VertexOwnInclusion())
        return 0;

    std::vector<vertex*> NpV = FindCloseVertices(pvertex);
    double dE = 0;
    // using modern C++ range-based loop with auto
    for (auto* v : NpV) {
            
            dE += CalculateNonbondedInteractionBetweenTwoVertices(v, pvertex);
    }
    
    return dE;  // Return computed nonbonded interaction energy
}

double InteractionBetweenInclusionsIn3D::CalculateNonbondedInteractionBetweenTwoVertices(vertex* v1, vertex* v2) const {

    if (!v1->VertexOwnInclusion() || !v2->VertexOwnInclusion()) return 0;

    double E;
    
   /* double Dist2 = v1->SquareDistanceFromAVertex(v2);
    if (Dist2 == 0) return 0;  // Prevent division by zero
    if (m_R0_2 < Dist2) return 0;

    Vec3D Dirc1 = (v1->GetL2GTransferMatrix()) * (v1->GetInclusion()->GetLDirection());
    Vec3D Dirc2 = (v2->GetL2GTransferMatrix()) * (v2->GetInclusion()->GetLDirection());

    double angle = Vec3D::dot(Dirc1, Dirc2);
    double E = -m_EP * (angle * angle - m_DAngle) / Dist2;*/

    return E;  // Return computed energy
}

std::vector<vertex*> InteractionBetweenInclusionsIn3D::FindCloseVertices(vertex* pV) const {
    std::vector<vertex*> pNV;  // Vector to store nearby vertices (excluding direct neighbors)
    
    // Get directly connected neighbor vertices
    std::vector<vertex*> npvertex = pV->GetVNeighbourVertex();

    // Use an unordered_set for efficient lookup when filtering
    std::unordered_set<vertex*> npSet(npvertex.begin(), npvertex.end());

    // Loop through all 26 neighboring voxel cells around the vertex's current voxel
    for (int n = -1; n < 2; n++) {
        for (int m = -1; m < 2; m++) {
            for (int s = -1; s < 2; s++) {
                
                // Get all vertices inside this neighboring cell
                std::vector<vertex*> CV = pV->GetVoxel()->GetANeighbourCell(n, m, s)->GetContentObjects();

                // Add to pNV only if the vertex is NOT in the direct neighbors set
                for (vertex* v : CV) {
                    if (npSet.find(v) == npSet.end()) {
                        pNV.push_back(v);
                    }
                }
            }
        }
    }

    return pNV;  // Return the filtered list of nearby, non-directly-connected vertices
}

// Implementation of CurrentState function
std::string InteractionBetweenInclusionsIn3D::CurrentState() const {
    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    state += " " + m_Input_Data;
    return state;
}



