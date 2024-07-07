#if !defined(AFX_CouplingtoFixedGlobalCurvature_H_334B21B8_C13C_2248_BF23_124095086255__INCLUDED_)
#define AFX_CouplingtoFixedGlobalCurvature_H_334B21B8_C13C_2248_BF23_124095086255__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This object affects the simulation results in every step. It is a global coupling, penalize the deviations of total mean curvature  from a target value
 It can also model the energy associated with the differnce in the inner and outer monolayer area
 The energy coupling is as below
 E=k/(2A)*(sum(2h-gc0)av)^2
 E=k/(2Ah^2)*(DA-DA0)^2
 */
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
class CouplingtoFixedGlobalCurvature
{
public:
    CouplingtoFixedGlobalCurvature();
    CouplingtoFixedGlobalCurvature(bool state);
    CouplingtoFixedGlobalCurvature(bool state,  double K, double GlobalC0);
    ~CouplingtoFixedGlobalCurvature();

       inline bool GetState()                           {return m_State;} // if the coupling is active
       inline double GetEnergy()                           {return m_Energy;} // Only part of energy asscoiated with this class
    
public:
    void CalculateGlobalVariables(std::vector<vertex *> Ver);
    double CalculateEnergyChange(double DA, double DC);
    void UpdateEnergyChange(double DA, double DC);

private:
    double m_K;        // energy coupling constant (in the constructor it will be devided by 2)
    double m_gC0;    // Global curvature
    bool m_State;   // to check if the coupling is active
    double m_TotalArea;     // total membrane area
    double m_TotalRCurvature;    //   SUM 2H_v*A_v
    double m_Energy;            // energy
};


#endif
