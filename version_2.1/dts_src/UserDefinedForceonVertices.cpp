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
    m_pV1 = m_pallV[Nfunction::String_to_Int(data[0])];
    m_pV2 = m_pallV[Nfunction::String_to_Int(data[1])];
    m_pV3 = m_pallV[Nfunction::String_to_Int(data[2])];
    m_K = Nfunction::String_to_Double(data[3]);
}
double UserDefinedForceonVertices::Energy_of_Force(vertex *pv, Vec3D dx) {
    

    Vec3D Force = CalculateForce(pv);
    double E = -Vec3D::dot(dx , Force);
    
    return E;
}
Vec3D UserDefinedForceonVertices::Force(vertex *pv) {
    Vec3D force = CalculateForce(pv);
    return force;
}
Vec3D UserDefinedForceonVertices::CalculateForce(vertex *pv) { // gives force in the local coordinate

    if(m_pV1 != pv && m_pV2 != pv && m_pV3 != pv){
        return 0;
    }

    Vec3D X0 = pv->GetPos();
    Vec3D X1 = m_pV1->GetPos();
    Vec3D X2 = m_pV2->GetPos();
    Vec3D X3 = m_pV3->GetPos();
    Vec3D f1 = X0 - X1;
    Vec3D f2 = X0 - X2;
    Vec3D f3 = X0 - X3;
    f1.normalize();
    f2.normalize();
    f3.normalize();
    Vec3D f = f1 + f2 + f3;
    f = f*m_K;
    return f;
}
std::string UserDefinedForceonVertices::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + m_Inputs;
    return state;
}

