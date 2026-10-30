#if !defined(AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Coupling the system energy to a potential for having osmotic pressure.
*/

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractVolumeCoupling.h"

class State;

class Apply_Osmotic_Pressure : public AbstractVolumeCoupling {
public:
    // Constructor
    Apply_Osmotic_Pressure(VAHGlobalMeshProperties *VHA, double gamma, double P0);

    // Destructor
    ~Apply_Osmotic_Pressure() override;

    // Initialize the volume coupling by calculating initial volumes and areas
    void Initialize(State* pstate) override;

    // Retrieve the current state of the volume coupling as a string
    std::string CurrentState() override;

    // Return the derived class name for identification purposes
    std::string GetDerivedDefaultReadName() const override {
        return "OsmoticPressure";
    }

    // Return the default read name for registration / identification
    inline static std::string GetDefaultReadName() {
        return "OsmoticPressure";
    }

    // Calculate the energy contribution from the volume coupling
    double GetCouplingEnergy() override;

    // Calculate the change in energy due to changes in area and volume
    double GetEnergyChange(double oldarea, double oldvolume,
                           double newarea, double newvolume) override;

private:
    // Compute the energy contribution from a given volume and area
    double Energy(double volume, double area);

    // Pointer to the state class
    State *m_pState;

    // Reference pressure
    double m_P0;

    // Initial volume
    double m_V0;

    // Osmotic pressure coefficient
    double m_Gamma;

    // Precomputed constant 1 / (6 * sqrt(pi))
    double m_6SQPI;
};

#endif // AFX_Apply_Osmotic_Pressure_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
