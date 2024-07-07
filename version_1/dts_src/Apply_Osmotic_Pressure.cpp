


#include "Apply_Osmotic_Pressure.h"
#include "Nfunction.h"
Apply_Osmotic_Pressure::Apply_Osmotic_Pressure()
{
    m_State = false;

}
Apply_Osmotic_Pressure::Apply_Osmotic_Pressure(bool State, std::string type, int eqsteps, double gamma,  double P0)
{
    double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2    
    m_TotalVolume = 0;
    m_TotalArea = 0;
    m_NoEQStep = eqsteps;
    m_P0 = P0;
    m_State = State;
    m_Gamma  = gamma;
    
    

    
}

Apply_Osmotic_Pressure::~Apply_Osmotic_Pressure()
{
    
}


double Apply_Osmotic_Pressure::VolumeofTrianglesAroundVertex(vertex *pv)
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
void Apply_Osmotic_Pressure::Initialize(std::vector<triangle *> pTriangle)
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
   double A0 = double(pTriangle.size())*3.0*sqrt(3)/4.0; /// a_t = sqrt(3)/4.0*l^2 ;; l=sqrt(1-3)

   m_V0 = m_6SQPI*sqrt(A0*A0*A0);
   
   //   std::cout<<"V0 "<<m_V0<<"  "<<V<< "\n";
     //       std::cout<<"V0 "<<A0<<"  "<<A<< "\n";
            


}
double Apply_Osmotic_Pressure::SingleTriangleVolume(triangle *pt)
{
    double vol=0;
    
    vertex* pv= pt->GetV1();
    double area= pt->GetArea();
    Vec3D Norm=pt->GetNormalVector();
    
    Vec3D R (pv->GetVXPos(),pv->GetVYPos(),pv->GetVZPos());
    
    // If the system is PBC broken, then the volume will be wrong .....
    {
        vertex* pv2= pt->GetV2();
        Vec3D R2 (pv2->GetVXPos(),pv2->GetVYPos(),pv2->GetVZPos());
        R2=R2-R;
        vertex* pv3= pt->GetV3();
        Vec3D R3 (pv3->GetVXPos(),pv3->GetVYPos(),pv3->GetVZPos());
        R3=R3-R;
        
        if(R2.dot(R2,R2)>4 || R3.dot(R3,R3)>4)
        {
            Nfunction f;
            std::string sms="---> Error, the system crossed the PBC while using volume coupling; use a large box .. ";
            f.Write_One_LogMessage(sms);
            std::cout<<sms<<std::endl;
            exit(0);
        }
    }
    
    vol=area*(R.dot(Norm,R))/3.0;
    
    
    return vol;
}

double Apply_Osmotic_Pressure::GetEnergyChange(int step, double oldvolume, double newvolume)
{

    double alpha=1;
    if(step<m_NoEQStep)
        alpha= double(step)/double(m_NoEQStep);
        

	double dv =  newvolume - oldvolume;
       double DE = -alpha*m_P0*(m_V0*log(1+dv/m_TotalVolume)-dv);
         //    double DE = -alpha*m_P0*dv;
      // std::cout<<DE<<"\n";

    return DE;
}
void Apply_Osmotic_Pressure::UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume)
{
    m_TotalVolume+=newvolume-oldvolume;
    m_TotalArea+=newarea-oldarea;
}

