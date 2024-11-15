#include "UserDefinedForceonVertices.h"
#include "State.h"

UserDefinedForceonVertices::UserDefinedForceonVertices(State *pState, std::string data_stream) : m_Inputs(data_stream), m_pState(pState)
{



}
UserDefinedForceonVertices::~UserDefinedForceonVertices() {
    
}
void UserDefinedForceonVertices::Initialize(){
    
    std::vector<std::string> data = Nfunction::Split(m_Inputs);
    m_pallV = m_pState->GetMesh()->GetActiveV();
}
double UserDefinedForceonVertices::Energy_of_Force(vertex *pv, Vec3D dx) {
    

    Vec3D Force = CalculateForce(pv);
    double E = -Vec3D::dot(dx , Force);
    
    return E;
}
Vec3D UserDefinedForceonVertices::CalculateForce(vertex *pv) { // gives force in the local coordinate

    Vec3D f;
    return f;
}
std::string UserDefinedForceonVertices::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + m_Inputs;
    return state;
}

