#include "SphericalVertexSubstrate.h"


SphericalVertexSubstrate::SphericalVertexSubstrate(std::string datastream) {

    std::vector<std::string> input_data = Nfunction::split(datastream);
    
    if(input_data.size() < 5){
        std::cout<<"---> error: SphericalVertexSubstrate the input is not enough \n";
        exit(0);
    }
    m_AdhesionStrength = Nfunction::String_to_Double(input_data[0]);
    m_Radius = Nfunction::String_to_Double(input_data[1]);
    m_Center(0) = Nfunction::String_to_Double(input_data[2]);
    m_Center(1) = Nfunction::String_to_Double(input_data[3]);
    m_Center(2) = Nfunction::String_to_Double(input_data[4]);

}
SphericalVertexSubstrate::~SphericalVertexSubstrate() {
    
}
double SphericalVertexSubstrate::GetCouplingEnergy(vertex *pvertex) {

    Vec3D R = pvertex->GetPos() - m_Center;
    double d_R = R.norm()-m_Radius;
    double E = -m_AdhesionStrength * 1/(d_R*d_R+1);
    return E;
}
std::string SphericalVertexSubstrate::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_AdhesionStrength) +" "+ Nfunction::D2S(m_Radius)+" "+ Nfunction::D2S(m_Center(0))+" "+ Nfunction::D2S(m_Center(1))+" "+ Nfunction::D2S(m_Center(2));
    return state;
}
