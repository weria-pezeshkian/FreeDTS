#if !defined(AFX_AbstractVolumeCoupling_H)
#define AFX_AbstractVolumeCoupling_H

#include <iostream>
#include "VAHGlobalMeshProperties.h"

/*
=======================================================
Developed 2024 by Weria Pezeshkian
Weria Pezeshkian (weria.pezeshkian@gmail.com)
Copyright (c) Weria Pezeshkian

This abstract base class defines the interface for volume coupling mechanisms
in a simulation environment. The primary purpose of this class is to couple
the system's energy to a volume constraint, allowing for dynamic adjustments
based on predefined parameters. Derived classes must implement the specific
details of these coupling mechanisms.
=======================================================
*/
class AbstractVolumeCoupling  {
public:
    // Constructor that initializes the base class with global mesh properties
    AbstractVolumeCoupling(VAHGlobalMeshProperties *pVHA) : m_pVAH(pVHA), m_TotalVolume(pVHA->m_TotalVolume), m_TotalArea(pVHA->m_TotalArea) {
    }

    // Virtual destructor to ensure proper cleanup of derived class objects
    virtual ~AbstractVolumeCoupling() {
    }

    // Pure virtual methods that must be implemented by derived classes

    // Initialize the volume coupling parameters and state
    virtual void Initialize(State* pState) = 0;

    // Retrieve the current state of the volume coupling as a string
    virtual std::string CurrentState() = 0;

    // Calculate the volume and area associated with a vertex ring
   // virtual void CalculateVolumeOfAVertexRing(vertex *pVertex, double &vol, double &area) = 0;

    // Calculate the volume and area associated with link triangles
  //  virtual void CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area) = 0;

    // Calculate the energy contribution from the volume coupling
    virtual double GetCouplingEnergy() = 0;

    // Calculate the change in energy due to changes in area and volume
    virtual double GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume) = 0;

    // Update the total area and volume based on changes
    virtual void UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume) = 0;

    // Retrieve the default name of the derived class for logging or identification purposes
    virtual inline std::string GetDerivedDefaultReadName() = 0;

    // Retrieve the base class default name for logging or identification purposes
    inline static std::string GetBaseDefaultReadName() { return "VolumeCoupling"; }

protected:
    VAHGlobalMeshProperties *m_pVAH;
    
    double &m_TotalVolume;
    double &m_TotalArea;
};

// Class representing no volume coupling, used as a default or placeholder
class NoCoupling : public AbstractVolumeCoupling {
public:
    // Constructor that initializes the base class with global mesh properties
    NoCoupling(VAHGlobalMeshProperties *VHA) : AbstractVolumeCoupling(VHA) {
    }

    // Destructor
    ~NoCoupling() {
    }

    // Return the derived class name for identification purposes
    virtual inline std::string GetDerivedDefaultReadName() { return "No"; }

    // Implementations of the virtual methods with no operation (no coupling)
    void Initialize(State* pState) { return; }
  //  void CalculateVolumeOfAVertexRing(vertex *pVertex, double &vol, double &area) { return; }
  //  void CalculateVolumeofALinkTriangles(links *p_link, double &vol, double &area) { return; }
    double GetCouplingEnergy() { return 0; }
    double GetEnergyChange(double oldarea, double oldvolume, double newarea, double newvolume) { return 0; }
    void UpdateArea_Volume(double oldarea, double oldvolume, double newarea, double newvolume) { return; }

    // Retrieve the current state as a string
    std::string CurrentState() {
        return GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    }
};

#endif
