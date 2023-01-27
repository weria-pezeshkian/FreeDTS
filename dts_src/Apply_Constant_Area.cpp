#include "Apply_Constant_Area.h"
#include "Nfunction.h"


Apply_Constant_Area::Apply_Constant_Area()
{
    m_State = false;

}
Apply_Constant_Area::Apply_Constant_Area(bool State, int eqsteps, double gamma,  double K0)
{
    double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2    
    m_TotalArea = 0;
    m_NoEQStep = eqsteps;
    m_K0 = K0;
    m_State = State;
    m_Gamma  = gamma;
    m_K0 = m_K0/2;
    m_NT = 1;
    
    if(m_Gamma<0 || m_Gamma>1)
    {
    std::cout<<"---> error in constant area; gamma is bad; make sure you know what are you doing \n";
    exit(0);
    }

    
}
Apply_Constant_Area::~Apply_Constant_Area()
{
    
}
double Apply_Constant_Area::AreaofTrianglesAroundVertex(vertex *pv)
{
    double A=0.0;
    
    
    std::vector <triangle *> pvT=pv->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = pvT.begin() ; it != pvT.end(); ++it)
		A+=(*it)->GetArea();

    
    return A;
}
//==========================================================
void Apply_Constant_Area::Initialize(std::vector<triangle *> pTriangle)
{
    double A=0.0;
	m_NT = double(pTriangle.size());
    for (std::vector<triangle *>::iterator it = pTriangle.begin() ; it != pTriangle.end(); ++it)
        A+=(*it)->GetArea();

    m_TotalArea = A;

   double A0 = double(pTriangle.size())*sqrt(3)/4.0; /// a_t = sqrt(3)/4.0*l^2 ;; l=sqrt(1-3)
   
   
   m_A0 = (1+2*m_Gamma)*A0;   // selecting between 1-3
    
    m_K0 = m_K0/m_NT;       // making the K0 t dependent
   
}
double Apply_Constant_Area::GetEnergyChange(int step, double oldarea, double newarea)
{

    double alpha=1;
    if(step<m_NoEQStep)
        alpha= double(step)/double(m_NoEQStep);
        
        double da = newarea - oldarea;
        
        //DE = (A+DA-A0)^2-(A-A0)^2 = 2(A-A0+DA/2)DA


       double DE = 2*(m_TotalArea+da/2-m_A0)*da;
       
       DE = alpha*m_K0*DE;



    return DE;
}
void Apply_Constant_Area::UpdateArea(double oldarea,  double newarea)
{

    m_TotalArea+=newarea-oldarea;
}


