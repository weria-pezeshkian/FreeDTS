
#include "VolumeCouplingSecondOrder.h"
#include "State.h"
#include "Nfunction.h"

// Constructor
VolumeCouplingSecondOrder::VolumeCouplingSecondOrder(VAHGlobalMeshProperties *VHA,  double DeltaP,  double K, double targetV)
    : AbstractVolumeCoupling(VHA)  {
    m_KV = K / 2;  // Halve the input parameter K
    m_TargetV = targetV;  // Set the target volume
    m_DeltaP = DeltaP;  // Set the pressure difference
    double pi = acos(-1);  // Calculate Pi
    m_6SQPI = 1.0 / (6.0 * sqrt(pi));  // Calculate 1 / (6 * sqrt(Pi))
        
}

VolumeCouplingSecondOrder::~VolumeCouplingSecondOrder() {
    
}
// Initialize the volume and total area
void VolumeCouplingSecondOrder::Initialize(State* pstate){

    m_pState = pstate;
    m_TotalVolume = 0;
    m_TotalArea = 0;
    std::vector<triangle *> all_tri = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
        m_TotalVolume += m_pVAH->CalculateSingleTriangleVolume(*it);
        m_TotalArea += (*it)->GetArea();
    }

    return;
}/*
void VolumeCouplingSecondOrder::CalculateVolumeOfAVertexRing(vertex *p_vertex, double &vol, double &area) {
    vol = 0.0;
    area = 0.0;
    const std::vector<triangle *>& ring_triangles = p_vertex->GetVTraingleList();
    for (std::vector<triangle *>::const_iterator it = ring_triangles.begin() ; it != ring_triangles.end(); ++it){
        vol+=CalculateSingleTriangleVolume(*it);
        area += (*it)->GetArea();

    }
    return;
}
void VolumeCouplingSecondOrder::CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area) {
    vol = 0.0;
    area = 0.0;
    
    vol += CalculateSingleTriangleVolume(p_link->GetTriangle());
    area += p_link->GetTriangle()->GetArea();
    
    vol += CalculateSingleTriangleVolume(p_link->GetMirrorLink()->GetTriangle());
    area += p_link->GetMirrorLink()->GetTriangle()->GetArea();

    return;
}

//==========================================================

double VolumeCouplingSecondOrder::CalculateSingleTriangleVolume(triangle *pTriangle){
    
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
double VolumeCouplingSecondOrder::GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume){

    double E1 =Energy(m_TotalVolume,m_TotalArea);
    double E2 =Energy(m_TotalVolume+newvolume-oldvolume,m_TotalArea+newarea-oldarea);
    return E2-E1;
}
double VolumeCouplingSecondOrder::Energy(double volume, double area){
   // m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2
        double E=0;
        double SQA = sqrt(area);
        double v0 = m_6SQPI*SQA*SQA*SQA;
        double v =volume/v0;
        return  -m_DeltaP*volume+m_KV*(v-m_TargetV)*(v-m_TargetV);
}
void VolumeCouplingSecondOrder::UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume) {
    m_TotalVolume += newvolume-oldvolume;
    m_TotalArea += newarea-oldarea;
}
std::string VolumeCouplingSecondOrder::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state +=  " " + Nfunction::D2S(m_DeltaP) + " " + Nfunction::D2S(2*m_KV) + " " + Nfunction::D2S(m_TargetV) ;
    return state;
}
double VolumeCouplingSecondOrder::GetCouplingEnergy(){
    return Energy(m_TotalVolume, m_TotalArea);
}
