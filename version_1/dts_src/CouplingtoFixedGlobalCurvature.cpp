


#include "CouplingtoFixedGlobalCurvature.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object affects the simulation results in every step. It is a global coupling, penalize the deviations of total mean curvature  from a target value
 It can also model the energy associated with the differnce in the inner and outer monolayer area
 The energy coupling is as below
 E=k/(2A)*(sum(2h-gc0)av)^2
 E=k/(2Ah^2)*(DA-DA0)^2
 */
CouplingtoFixedGlobalCurvature::CouplingtoFixedGlobalCurvature()
{
    m_State = false;
    m_TotalArea = 0;
    m_TotalRCurvature = 0;
    m_Energy = 0;
    m_K = 0;
}
CouplingtoFixedGlobalCurvature::CouplingtoFixedGlobalCurvature(bool state)
{
    m_K = 0;
    m_gC0 = 0;
    m_State = state;
    m_TotalArea = 0;
    m_TotalRCurvature = 0;
    m_Energy = 0;

}
CouplingtoFixedGlobalCurvature::CouplingtoFixedGlobalCurvature(bool state,  double Gkappa, double GlobalC0)
{
    m_K = Gkappa/2.0;
    m_gC0 = GlobalC0;
    m_State = state;
    m_TotalArea = 0;
    m_TotalRCurvature = 0;
    m_Energy = 0;
}

CouplingtoFixedGlobalCurvature::~CouplingtoFixedGlobalCurvature()
{
    
}
// A function to initialize total area, total curvature and energy asscoiated with the coupling
void CouplingtoFixedGlobalCurvature::CalculateGlobalVariables(std::vector<vertex *> Ver)
{
    double ta = 0;
    double tc=0;
    for (std::vector<vertex *>::iterator it = Ver.begin() ; it != Ver.end(); ++it)
    {
        std::vector<double> Curve=(*it)->GetCurvature();
        double a = (*it)->GetArea();
        ta+=a;
        tc+=(Curve.at(0)+Curve.at(1))*a;
    }
   m_TotalArea = ta;
   m_TotalRCurvature = tc;
    
    double dh=(m_TotalRCurvature-m_gC0*m_TotalArea);
    m_Energy =m_K/(m_TotalArea)*dh*dh;
}
// when a vertex moves or a link flips or any other changes, we can see the changes in the total area and total mean curavture:
// this function can tell us how much this changes cost energy but it does not update the change since it can be rejected.
double CouplingtoFixedGlobalCurvature::CalculateEnergyChange(double DA, double DC)
{
    double de=0;
    double dh=(m_TotalRCurvature+DC-m_gC0*(m_TotalArea+DA));
    double e= m_K/((m_TotalArea+DA))*dh*dh;  // new energy
    
    de=e-m_Energy;   // de = newenergy -oldenergy
    
    return de;
}
// this function update the changes if the move get accetped
void CouplingtoFixedGlobalCurvature::UpdateEnergyChange(double DA, double DC)
{
    double de=0;
    
    m_TotalArea+= DA;
    m_TotalRCurvature+= DC;
    double dh=(m_TotalRCurvature-m_gC0*m_TotalArea);
    m_Energy= m_K*dh*dh/(m_TotalArea);
 
}
