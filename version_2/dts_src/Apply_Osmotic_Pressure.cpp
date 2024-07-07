#include "Apply_Osmotic_Pressure.h"
#include "State.h"
#include "Nfunction.h"
/*
 
 Delta E_osmos = -RT[c_in*V_ini*Log (V/V_ini) - c_out * (V-V_ini)]
 
 Delta E_osmos = -m_P0*(m_V0*log(1+dv/m_TotalVolume)-dv);
 
 m_V0 = m_6SQPI*sqrt(A0*A0*A0);
 A0 = no_T*(minT area)*(1+2*gamma)^2

 
 */

Apply_Osmotic_Pressure::Apply_Osmotic_Pressure(VAHGlobalMeshProperties *VHA,  double gamma,  double P0) :
 AbstractVolumeCoupling(VHA) {
     
    m_P0 = P0;
    m_Gamma  = gamma;
     
     double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));
     
     if(m_Gamma<0 || m_Gamma>1) {
         std::cout<<"---> error for Osmotic Pressure; gamma is bad; make sure you know what are you doing \n";
         exit(0);
     }

}

Apply_Osmotic_Pressure::~Apply_Osmotic_Pressure()
{
    
}
void Apply_Osmotic_Pressure::Initialize(State* pstate){
    m_pState = pstate;
    m_TotalVolume = 0;
    m_TotalArea = 0;
    m_CalculatedGlobalVariable = true;
    std::vector<triangle *> all_tri = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
        m_TotalVolume += m_pVAH->CalculateSingleTriangleVolume(*it);
        m_TotalArea += (*it)->GetArea();
    }
    
    double l2 = 1+2*m_Gamma;
    double A0 = double(all_tri.size())*sqrt(3)/4.0*l2; /// a_t = sqrt(3)/4.0*l^2 ;; l=sqrt(1-3)

    m_V0 = m_6SQPI*sqrt(A0*A0*A0);

    return;
}/*
void Apply_Osmotic_Pressure::CalculateVolumeOfAVertexRing(vertex *p_vertex, double &vol, double &area) {
    vol = 0.0;
    area = 0.0;
    const std::vector<triangle *>& ring_triangles = p_vertex->GetVTraingleList();
    for (std::vector<triangle *>::const_iterator it = ring_triangles.begin() ; it != ring_triangles.end(); ++it){
        vol+=CalculateSingleTriangleVolume(*it);
        area += (*it)->GetArea();

    }
    return;
}
void Apply_Osmotic_Pressure::CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area) {
    vol = 0.0;
    area = 0.0;
    
    vol += CalculateSingleTriangleVolume(p_link->GetTriangle());
    area += p_link->GetTriangle()->GetArea();
    
    vol += CalculateSingleTriangleVolume(p_link->GetMirrorLink()->GetTriangle());
    area += p_link->GetMirrorLink()->GetTriangle()->GetArea();

    return;
}
double Apply_Osmotic_Pressure::CalculateSingleTriangleVolume(triangle *pTriangle){
    
    if(m_pState->GetMesh()->GetHasCrossedPBC()){
        *(m_pState->GetTimeSeriesLog()) << "---> the system has crossed the PBC while volume is being calculated.";
        *(m_pState->GetTimeSeriesLog()) << " SOLUTION: Restart the simulation and center the system. Also, activate the command for centering the box.";

         exit(0);
    }
    
    double T_area = pTriangle->GetArea();
    Vec3D Normal_v = pTriangle->GetNormalVector();
    Vec3D Pos = pTriangle->GetV1()->GetPos();

    // Compute triangle volume
    return T_area * (Pos.dot(Normal_v, Pos)) / 3.0;
}*/
double Apply_Osmotic_Pressure::GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume){

    double dv =  newvolume - oldvolume;
    return -m_P0*(m_V0*log(1+dv/m_TotalVolume)-dv);
}
void Apply_Osmotic_Pressure::UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume) {
    m_TotalVolume += newvolume-oldvolume;
    m_TotalArea += newarea-oldarea;
}
std::string Apply_Osmotic_Pressure::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " " + Nfunction::D2S(m_P0) + " " + Nfunction::D2S(m_Gamma);
    return state;
}
double Apply_Osmotic_Pressure::GetCouplingEnergy(){
    return Energy(m_TotalVolume, m_TotalArea);
}
double Apply_Osmotic_Pressure::Energy(double volume, double area){
    std::cout<<" this function has not been completed 234 \n";
    return 0;
}
