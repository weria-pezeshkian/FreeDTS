


#include "CmdVolumeCouplingSecondOrder.h"
CmdVolumeCouplingSecondOrder::CmdVolumeCouplingSecondOrder()
{
}
CmdVolumeCouplingSecondOrder::CmdVolumeCouplingSecondOrder(bool State, int eqsteps, double DeltaP,  double K, double targetV)
{
    double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2

    
    
    m_TotalVolume = 0;
    m_TotalArea = 0;
    m_NoEQStep = eqsteps;
    m_KV = K/2;
    m_State = State;
    m_TargetV = targetV;
    m_DeltaP  = DeltaP;
    
}

CmdVolumeCouplingSecondOrder::~CmdVolumeCouplingSecondOrder()
{
    
}


double CmdVolumeCouplingSecondOrder::VolumeofTrianglesAroundVertex(vertex *pv)
{
    double vol=0;
    
    
    std::vector <triangle *> pvT=pv->GetVTraingleList();

    
    for (std::vector<triangle *>::iterator it = pvT.begin() ; it != pvT.end(); ++it)
    {
        
        vol+=SingleTriangleVolume(*it);
        
    }
    
    return vol;
}
//==========================================================
void CmdVolumeCouplingSecondOrder::Initialize(std::vector<triangle *> pTriangle)
{
    double V=0.0;
    double A=0.0;

    for (std::vector<triangle *>::iterator it = pTriangle.begin() ; it != pTriangle.end(); ++it)
    {
        V+=SingleTriangleVolume(*it);
        A+=(*it)->GetArea();
    }
    m_TotalVolume = V;
    m_TotalArea = A;
}
double CmdVolumeCouplingSecondOrder::SingleTriangleVolume(triangle *pt)
{
    double vol=0;
    
    vertex* pv= pt->GetV1();
    double area= pt->GetArea();
    Vec3D Norm=pt->GetNormalVector();
    
    Vec3D R (pv->GetVXPos(),pv->GetVYPos(),pv->GetVZPos());
    
    
    
    vol=area*(R.dot(Norm,R))/3.0;
    
    
    return vol;
}

double CmdVolumeCouplingSecondOrder::GetEnergyChange(int step, double oldarea, double oldvolume, double newarea, double newvolume)
{
    double DE=0;
    double alpha=m_TargetV;
    if(step<m_NoEQStep)
        alpha= 1-(1-alpha)*double(step)/double(m_NoEQStep);
        

        double E1 =Energy(m_TotalVolume,m_TotalArea,alpha);
        double E2 =Energy(m_TotalVolume+newvolume-oldvolume,m_TotalArea+newarea-oldarea,alpha);
        DE = E2-E1;
    

    return DE;
}
double CmdVolumeCouplingSecondOrder::Energy(double volume, double area, double alpha)
{
    
    
        double E=0;
        double SQA = sqrt(area);
        double v0 = m_6SQPI*SQA*SQA*SQA;
        double v =volume/v0;
        E = -m_DeltaP*volume+m_KV*(v-alpha)*(v-alpha)*v0*v0;
        
        return E;
}
void CmdVolumeCouplingSecondOrder::UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume)
{
    m_TotalVolume+=newvolume-oldvolume;
    m_TotalArea+=newarea-oldarea;
}

