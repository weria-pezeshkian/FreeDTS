#if !defined(AFX_Apply_Constant_Vertex_Area_H_994B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Apply_Constant_Vertex_Area_H_994B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class is created does not work yet.
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
class Apply_Constant_Vertex_Area
{
public:
    Apply_Constant_Vertex_Area();
    Apply_Constant_Vertex_Area(bool State,int eqsteps, double gamma,  double P0);
    ~Apply_Constant_Vertex_Area();



       inline double GetTotalArea()                  {return m_TotalArea;}
       inline bool GetState()                   {return m_State;}
    


public:
    

    
    
    
    double AreaofTrianglesAroundVertex(vertex * pVeretx);   ///
    void Initialize(std::vector<triangle *> pTriangle);   ///
    double GetEnergyChange(int step, double oldarea,  double newarea);
    void UpdateArea(double oldarea, double newarea);

    //=====
    
    
    
private:
    double m_TotalArea;
    int m_NoEQStep;
    double m_K0;
    double m_A0;
    double m_NT;  // number of the triangles in the system
    bool m_State;
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
