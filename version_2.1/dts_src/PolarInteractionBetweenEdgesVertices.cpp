#include "PolarInteractionBetweenEdgesVertices.h"
#include "State.h"

// Constructor
PolarInteractionBetweenEdgesVertices::PolarInteractionBetweenEdgesVertices(State* pstate, std::string input_data) {
    m_Input_Data = input_data;
    m_pState = pstate;
}
PolarInteractionBetweenEdgesVertices::~PolarInteractionBetweenEdgesVertices() {
    // Destructor
}
void PolarInteractionBetweenEdgesVertices::Initialize() {
    std::vector<std::string> data = Nfunction::Split(m_Input_Data);
    
    // Prevent out-of-bounds access
    if (data.size() < 3) {
        std::cout << "---> error: insufficient input data for "
                  << this->GetDefaultReadName() << ".\n";
        exit(-1);
    }

    // Use correct parsing function (double instead of int)
    m_EP = Nfunction::String_to_Double(data[0]);
    m_R0 = Nfunction::String_to_Double(data[1]);
    m_DAngle = Nfunction::String_to_Double(data[2]);

    m_pEdgeV = m_pState->GetMesh()->GetEdgeV();
    m_pBox = m_pState->GetMesh()->GetBox();
    m_R0_2 = m_R0 * m_R0;

    // Improved error messages and checks
    if (m_EP < 0) {
        std::cout << "---> error: assigned interaction strength in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }
    if (m_DAngle < 0 || m_DAngle > 90) {
        std::cout << "---> error: assigned angle in "
                  << this->GetDefaultReadName() << " should be between 0 and 90 degrees.\n";
        exit(-1);
    }
    if (m_R0 < 0) {
        std::cout << "---> error: assigned radius in "
                  << this->GetDefaultReadName() << " should be greater than zero.\n";
        exit(-1);
    }

    // More precise and readable angle conversion
    m_DAngle = cos(m_DAngle * (M_PI / 180.0));
    m_DAngle *= m_DAngle;  // Squaring it

    return;
}
double PolarInteractionBetweenEdgesVertices::GetVertexNonBondedEnergy(vertex* pvertex) const {
    
    if (pvertex->GetVertexType() != 1)
        return 0;
    
    // Check if the vertex owns an inclusion
    if (!pvertex->VertexOwnInclusion())
        return 0;

    // Prevent potential crashes if m_pEdgeV is empty
    if (m_pEdgeV.empty()) return 0;

    double dE = 0;
    
    // using modern C++ range-based loop with auto
    for (auto* v : m_pEdgeV) {
        if (v != pvertex &&
            v->GetEdgeLink() && v->GetEdgeLink()->GetV2() != pvertex &&
            v->GetPrecedingEdgeLink() && v->GetPrecedingEdgeLink()->GetV1() != pvertex) {
            
            dE += CalculateNonbondedInteractionBetweenTwoVertices(v, pvertex);
        }
    }
    
    return dE;  // Return computed nonbonded interaction energy
}

double PolarInteractionBetweenEdgesVertices::CalculateNonbondedInteractionBetweenTwoVertices(vertex* v1, vertex* v2) const {

    if (!v1->VertexOwnInclusion() || !v2->VertexOwnInclusion()) return 0;

    double Dist2 = v1->SquareDistanceFromAVertex(v2);
    if (Dist2 == 0) return 0;  // Prevent division by zero
    if (m_R0_2 < Dist2) return 0;

    Vec3D Dirc1 = (v1->GetL2GTransferMatrix()) * (v1->GetInclusion()->GetLDirection());
    Vec3D Dirc2 = (v2->GetL2GTransferMatrix()) * (v2->GetInclusion()->GetLDirection());

    double angle = Vec3D::dot(Dirc1, Dirc2);
    double E = -m_EP * (angle * angle - m_DAngle) / Dist2;

    return E;  // Return computed energy
}
// Implementation of CurrentState function
std::string PolarInteractionBetweenEdgesVertices::CurrentState() const {
    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    state += " " + m_Input_Data;
    return state;
}



