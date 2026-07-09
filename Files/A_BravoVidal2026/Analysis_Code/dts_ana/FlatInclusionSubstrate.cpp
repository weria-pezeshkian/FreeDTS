#include "FlatInclusionSubstrate.h"


FlatInclusionSubstrate::FlatInclusionSubstrate(std::string datastream) {

    std::vector<std::string> input_data = Nfunction::split(datastream);
    
    if(input_data.size() < 3){
        std::cout<<"---> error: FlatInclusionSubstrate the input is not enough \n";
        exit(0);
    }
    m_AdhesionStrength = Nfunction::String_to_Double(input_data[0]);
    m_Z = Nfunction::String_to_Double(input_data[1]);
    m_inctype = input_data[2];
    std::cout<<m_inctype<<std::endl;


}
FlatInclusionSubstrate::~FlatInclusionSubstrate() {
    
}
double FlatInclusionSubstrate::GetCouplingEnergy(vertex *pvertex) {


    if(!pvertex->VertexOwnInclusion()) {

        //---> if the vertex is just a membrane; 
        double E=0;
        return E;
    }
    else{
        
        inclusion *p_inc=pvertex->GetInclusion();
        std::string inctype=p_inc->m_IncType->ITName;
        
        if (inctype==m_inctype){

            double d_Z = pvertex->GetZPos() - m_Z;
            double E = -m_AdhesionStrength * 1/(d_Z*d_Z+1);
            return E;

        }
        else {

            double E = 0;

            return E;}
    }
}
std::string FlatInclusionSubstrate::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_AdhesionStrength) +" "+ Nfunction::D2S(m_Z);
    return state;
}
