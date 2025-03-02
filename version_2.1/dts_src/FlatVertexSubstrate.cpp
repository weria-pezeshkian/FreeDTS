#include "FlatVertexSubstrate.h"


FlatVertexSubstrate::FlatVertexSubstrate(std::string datastream) {

    std::vector<std::string> input_data = Nfunction::split(datastream);
    
    if(input_data.size() < 3){
        std::cout<<"---> error: FlatVertexSubstrate the input is not enough \n";
        exit(0);
    }
    m_AdhesionStrength = Nfunction::String_to_Double(input_data[0]);
    m_Depth = Nfunction::String_to_Double(input_data[1]);
    m_Z = Nfunction::String_to_Double(input_data[2]);

    if(m_Depth <= 0){
        std::cout<<"---> error: FlatVertexSubstrate depth cannot be zero  \n";
        exit(0);
    }

}
FlatVertexSubstrate::~FlatVertexSubstrate() {
    
}
double FlatVertexSubstrate::GetCouplingEnergy(vertex *pvertex) {

    double d_Z = pvertex->GetZPos() - m_Z;
   // double E = -m_AdhesionStrength * 1/(d_Z*d_Z+1);
    double E = 0;
    if(fabs(d_Z)<3*m_Depth){
        d_Z *= d_Z;
        double r4 = d_Z*d_Z;
         E = m_AdhesionStrength * (r4-2*d_Z);
    }

    return E;
}
std::string FlatVertexSubstrate::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_AdhesionStrength) +" "+ Nfunction::D2S(m_Z);
    return state;
}
