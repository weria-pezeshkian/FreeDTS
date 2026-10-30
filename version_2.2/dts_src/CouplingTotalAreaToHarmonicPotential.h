#if !defined(AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 * CouplingTotalAreaToHarmonicPotential.h
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 *
 * Description:
 * This class implements a harmonic potential-based coupling that constrains
 * the total membrane surface area toward a target (reference) area.
 *
 * The energy contribution penalizes deviations of the system's total area
 * from the preferred value using a quadratic (harmonic) form:
 *
 *     E ~ (A - A0)^2
 *
 * where A is the instantaneous total area and A0 is the reference area.
 *
 * Functionality:
 * - Computes total area contributions from vertices and triangle links
 * - Evaluates energy changes due to area fluctuations
 * - Provides coupling energy for use in Monte Carlo / dynamic simulations
 *
 * Dependencies:
 * - VAHGlobalMeshProperties
 * - State
 * - vertex, triangle, links
 * - AbstractTotalAreaCoupling (base class)
 *
 * Notes:
 * Note, this is a global based potentials, and in MC models cannot be parallelized
 */

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractTotalAreaCoupling.h"


class State;
class CouplingTotalAreaToHarmonicPotential : public  AbstractTotalAreaCoupling {
public:
    CouplingTotalAreaToHarmonicPotential(VAHGlobalMeshProperties *VHA, double K0, double gamma  );
    ~CouplingTotalAreaToHarmonicPotential();


    void Initialize(State *pstate);   ///
    double CalculateEnergyChange(double oldarea,  double newarea);
    double GetCouplingEnergy();
    double CalculateAreaofALinkTriangles(links *p_link);
    double CalculateAreaOfAVertexRing(vertex * pVeretx);   ///
    inline  std::string GetDerivedDefaultReadName()  {return "HarmonicPotential";}
    inline static std::string GetDefaultReadName() {return "HarmonicPotential";}
    std::string CurrentState();

private:
    double m_K0;
    double m_A0;
    double m_NT;  // number of the triangles in the system
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

    State *m_pState;
};


#endif
