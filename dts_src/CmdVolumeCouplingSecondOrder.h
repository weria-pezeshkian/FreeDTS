#if !defined(AFX_CmdVolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_CmdVolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class is created in version 1.2 to couple the system to a force
for taraget CmdVolumeCouplingSecondOrder.
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
class CmdVolumeCouplingSecondOrder
{
public:
    CmdVolumeCouplingSecondOrder();
    CmdVolumeCouplingSecondOrder(bool State, int eqsteps, double DeltaP,  double K, double targetV);
    ~CmdVolumeCouplingSecondOrder();


       inline double GetTotalVolume()                  {return m_TotalVolume;}
       inline double GetTotalArea()                  {return m_TotalArea;}
       inline bool GetState()                   {return m_State;}
    


public:
    

    
    
    
    double VolumeofTrianglesAroundVertex(vertex * pVeretx);   /// this does not mean anything outside of this code
    double SingleTriangleVolume(triangle * ptriangle);   /// and this one
    void Initialize(std::vector<triangle *> pTriangle);   ///
    double GetEnergyChange(int step, double oldarea, double oldvolume, double newarea, double newvolume);
    double Energy(double volume, double area, double alpha);
    void UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume);

    //=====
    
    
    
private:
    double m_TotalVolume;
    double m_TotalArea;
    int m_NoEQStep;
    double m_KV;
    bool m_State;
    double m_TargetV;
    double m_DeltaP;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
