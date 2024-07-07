#include "CouplingTotalAreaToHarmonicPotential.h"
#include "Nfunction.h"
#include "State.h"


CouplingTotalAreaToHarmonicPotential::CouplingTotalAreaToHarmonicPotential(VAHGlobalMeshProperties *VHA,  double K0, double gamma) : AbstractTotalAreaCoupling(VHA) {
    
    m_Gamma  = gamma;
    m_K0 = K0/2;
    
    double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2    

    
    if(m_Gamma<0 || m_Gamma>1) {
        std::cout<<"---> error in constant area; gamma is bad; make sure you know what are you doing \n";
        exit(0);
    }
}
CouplingTotalAreaToHarmonicPotential::~CouplingTotalAreaToHarmonicPotential() {
    
}
void CouplingTotalAreaToHarmonicPotential::Initialize(State* pstate){
    m_pState = pstate;
    m_TotalArea = 0;
    m_CalculatedGlobalVariable = true;
    std::vector<triangle *> all_tri = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
        m_TotalArea += (*it)->GetArea();
    }
    
    m_A0 = (1+2*m_Gamma)*double(all_tri.size())*sqrt(3)/4.0;   // selecting between 1-3
    m_K0 = m_K0 * double(all_tri.size())/(m_A0 * m_A0); // E=0.5*N*K*(A/A0-1)^2; K=0.5*N*K/A0^2==> E=0.5*N*K/A0^2(A-A0)^2
    
    return;
}

double CouplingTotalAreaToHarmonicPotential::CalculateEnergyChange(double oldarea, double newarea) {

    double dA = newarea - oldarea;
    double DE = (2*(m_TotalArea-m_A0)+dA) * dA;     // DE = alpha * (2 * [A - A0] +dA) * dA, where alpha = k/(2*A0^2)*NT
       
    return m_K0*DE;
}
double CouplingTotalAreaToHarmonicPotential::GetCouplingEnergy(){
    
    double E = 0;
    
    E = m_K0 * (m_TotalArea - m_A0) * (m_TotalArea - m_A0);
    
    return E;
}
void CouplingTotalAreaToHarmonicPotential::UpdateTotalArea(double oldarea,  double newarea) {
    m_TotalArea += newarea-oldarea;
}
double CouplingTotalAreaToHarmonicPotential::CalculateAreaOfAVertexRing(vertex *p_vertex) {
    double A=0.0;
    
    std::vector <triangle *> pvT = p_vertex->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = pvT.begin() ; it != pvT.end(); ++it){
        A += (*it)->GetArea();
    }

    return A;
}
double CouplingTotalAreaToHarmonicPotential::CalculateAreaofALinkTriangles(links *p_link){
     
    double  area = p_link->GetTriangle()->GetArea();
    area += p_link->GetMirrorLink()->GetTriangle()->GetArea();

    return area;
}
std::string CouplingTotalAreaToHarmonicPotential::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " " + Nfunction::D2S(m_K0) + " " + Nfunction::D2S(m_Gamma);
    return state;
}

