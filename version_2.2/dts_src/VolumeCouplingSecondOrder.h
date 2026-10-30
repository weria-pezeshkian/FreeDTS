#if !defined(AFX_VolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_VolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_

/*
=======================================================
Developed by Weria Pezeshkian
Weria Pezeshkian (weria.pezeshkian@gmail.com)
Copyright (c) Weria Pezeshkian

This class, created in version 1.2, extends the AbstractVolumeCoupling base class.
It couples the energy of the system to a volume constraint using a second-order polynomial approach.

Key parameters:
- DeltaP: Pressure difference applied to the system.
- K: Coupling strength parameter.
- targetV: Target reduced volume.

The class computes volume contributions of vertices and link triangles,
and calculates total energy changes due to volume modifications.
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
    VolumeCouplingSecondOrder(VAHGlobalMeshProperties *VHA, double DeltaP, double K, double targetV);

    // Destructor
    ~VolumeCouplingSecondOrder() override;

    // Initialize the volume coupling by calculating initial volumes and areas
    void Initialize(State* pstate) override;

    // Retrieve the current state of the volume coupling as a string
    std::string CurrentState() override;

    // Return the derived class name for identification purposes
    std::string GetDerivedDefaultReadName() const override {
        return "SecondOrder";
    }

    // Return the default read name for registration / identification
    inline static std::string GetDefaultReadName() {
        return "SecondOrder";
    }

    // Calculate the energy contribution from the volume coupling
    double GetCouplingEnergy() override;

    // Calculate the change in energy due to changes in area and volume
    double GetEnergyChange(double oldarea, double oldvolume,
                           double newarea, double newvolume) override;

    friend class NonequilibriumCommands;

private:
    // Compute the energy contribution from a given volume and area
    double Energy(double volume, double area);

    // Pointer to the state class
    State *m_pState;

    // Coupling parameter (half of input K)
    double m_KV;

    // Target reduced volume
    double m_TargetV;

    // Pressure difference
    double m_DeltaP;

    // Precomputed constant 1/(6*sqrt(pi))
    double m_6SQPI;
};

#endif // AFX_VolumeCouplingSecondOrder_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
