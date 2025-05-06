#include "AnisotropicStretchingByNematicField.h"


AnisotropicStretchingByNematicField::AnisotropicStretchingByNematicField(double K_p, double K_n) : m_Kp(K_p), m_Kn(K_n){

}
AnisotropicStretchingByNematicField::~AnisotropicStretchingByNematicField() {
    
}
void AnisotropicStretchingByNematicField::Initialize(){
    return;
}
double AnisotropicStretchingByNematicField::Energy(vertex *pv) {
    
    if(m_Kn == 0 && m_Kp == 0 ){
        return 0;
    }
    if(!pv->VertexOwnInclusion() ||  m_Kn == 0 ){
        return 0;
    }
    double En = 0;
    
    
   // Vec3D nem1_d = (v1->GetInclusion())->GetLDirection();   // this is in v1



    return En;
}
std::string AnisotropicStretchingByNematicField::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_Kp) +" "+ Nfunction::D2S(m_Kn) ;
    return state;
}

