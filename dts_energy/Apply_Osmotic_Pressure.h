#if !defined(AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class is created in version 1.2 to couple the system to a force
for taraget Apply_Osmotic_Pressure.
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
class Apply_Osmotic_Pressure
{
public:
    Apply_Osmotic_Pressure();
    Apply_Osmotic_Pressure(bool State, std::string type, int eqsteps, double gamma,  double P0);
    ~Apply_Osmotic_Pressure();


       inline double GetTotalVolume()                  {return m_TotalVolume;}
       inline double GetTotalArea()                  {return m_TotalArea;}
       inline bool GetState()                   {return m_State;}
    


public:
    

    
    
    
    double VolumeofTrianglesAroundVertex(vertex * pVeretx);   ///
    double SingleTriangleVolume(triangle * ptriangle);   ///
    void Initialize(std::vector<triangle *> pTriangle);   ///
    double GetEnergyChange(int step, double oldvolume,  double newvolume);
    void UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume);

    //=====
    
    
    
private:
    double m_TotalVolume;
    double m_TotalArea;
    int m_NoEQStep;
    double m_P0;
    bool m_State;
    double m_V0;
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
