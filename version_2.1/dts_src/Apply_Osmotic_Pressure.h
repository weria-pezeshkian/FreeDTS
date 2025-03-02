#if !defined(AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 coupling the system energy to a potential for having osmotic pressure.
 
 
 
 
 
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractVolumeCoupling.h"

class State;

class Apply_Osmotic_Pressure : public AbstractVolumeCoupling {
public:
    Apply_Osmotic_Pressure(VAHGlobalMeshProperties *VHA,  double gamma,  double P0);
    ~Apply_Osmotic_Pressure();


    // Initialize the volume coupling by calculating initial volumes and areas
    void Initialize(State* pstate);
    
    // Retrieve the current state of the volume coupling as a string
    std::string CurrentState();
    
    // Return the derived class name for identification purposes
    inline std::string GetDerivedDefaultReadName() { return "OsmoticPressure"; }
    
    // Return the default read name for identification purposes
    inline static std::string GetDefaultReadName() { return "OsmoticPressure"; }
    
    // Calculate the volume and area associated with a vertex ring
  //  void CalculateVolumeOfAVertexRing(vertex *pVertex, double &vol, double &area);
    
    // Calculate the volume and area associated with link triangles
   // void CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area);

    // Calculate the energy contribution from the volume coupling
    double GetCouplingEnergy();
    
    // Calculate the change in energy due to changes in area and volume
    double GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume);
    
    //=====
    
    
    
private:
    // Calculate the volume of a single triangle element
  //  double CalculateSingleTriangleVolume(triangle *pTriangle);
    
    // Compute the energy contribution from a given volume and area
    double Energy(double volume, double area);
    
    State *m_pState;      // pointer to the state class
    double m_P0;          // reference pressure 
    double m_V0;
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
