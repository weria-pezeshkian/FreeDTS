#if !defined(AFX_VolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_VolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
=======================================================
Developed by Weria Pezeshkian
Weria Pezeshkian (weria.pezeshkian@gmail.com)
Copyright (c) Weria Pezeshkian

This class, created in version 1.2, extends the AbstractVolumeCoupling base class.
It couples the energy of the system to a volume constraint using a second-order polynomial approach.

The class enforces a volume constraint by introducing a potential that keeps the system's volume
close to a target volume. The coupling mechanism utilizes both first and second-order polynomial terms
to accurately control the volume changes.

Key parameters include:
- DeltaP: The pressure difference applied to the system.
- K: The coupling parameter that controls the strength of the volume constraint.
- targetV: The target reduced volume the system should aim to maintain.

The class also computes the volume contributions of individual elements such as vertices and link triangles,
and calculates the total energy change due to volume modifications.

========================================================
*/
#include "AbstractVolumeCoupling.h"
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"


class State;

class VolumeCouplingSecondOrder : public AbstractVolumeCoupling {
public:
    // Constructor initializing the coupling parameters and target volume
    VolumeCouplingSecondOrder(VAHGlobalMeshProperties *VHA,  double DeltaP, double K, double targetV);
    
    // Destructor
    ~VolumeCouplingSecondOrder();

    // Initialize the volume coupling by calculating initial volumes and areas
    void Initialize(State* pstate);
    
    // Retrieve the current state of the volume coupling as a string
    std::string CurrentState();
    
    // Return the derived class name for identification purposes
    inline std::string GetDerivedDefaultReadName() { return "SecondOrder"; }
    
    // Return the default read name for identification purposes
    inline static std::string GetDefaultReadName() { return "SecondOrder"; }
    
    // Calculate the volume and area associated with a vertex ring
   // void CalculateVolumeOfAVertexRing(vertex *pVertex, double &vol, double &area);
    
    // Calculate the volume and area associated with link triangles
   // void CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area);

    // Calculate the energy contribution from the volume coupling
    double GetCouplingEnergy();
    
    // Calculate the change in energy due to changes in area and volume
    double GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume);
    
private:
    // Calculate the volume of a single triangle element
  //  double CalculateSingleTriangleVolume(triangle *pTriangle);
    
    // Compute the energy contribution from a given volume and area
    double Energy(double volume, double area);

    State *m_pState;      // pointer to the state class
    double m_KV;          // Coupling parameter (half of the input parameter K)
    double m_TargetV;     // Target reduced volume
    double m_DeltaP;      // Pressure difference
    double m_6SQPI;       // Precomputed constant 1/(6*sqrt(pi))
    


};


#endif
