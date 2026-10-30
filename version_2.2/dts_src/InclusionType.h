#ifndef INCLUSIONTYPE_H
#define INCLUSIONTYPE_H

#include "SimDef.h"
/*
 * @brief InclusionType object represents a type of inclusion in the simulation.
 *
 * This class defines the properties and characteristics of an inclusion type,
 * including its name, ID, symmetry, curvature parameters, line tension, and rigidity.
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */

class InclusionType {
    
public:
    InclusionType();
    /*
      * @brief Constructor for InclusionType.
      * @param name The name of the inclusion type.
      * @param id The ID of the inclusion type.
      * @param n The in-plane symmetry of the inclusion type.
      * @param k Type Kappa parameter.
      * @param kg KappaG parameter.
      * @param k1 K_|| parameter.
      * @param k2 K_norm parameter.
      * @param c0 Curvature parameter.
      * @param c1 Direction curvature 1 parameter.
      * @param c2 Direction curvature 2 parameter.
      * @param lam Line tension parameter.
      * @param ekg Geodesic rigidity parameter.
      * @param ekn Normal curvature line rigidity parameter.
      * @param ecn Spontaneous normal curvature parameter.
      */
    InclusionType(std::string, int id, int n, double k, double kg, double k1, double k2, double c0, double c1, double c2, double lam, double ekg, double ekn, double ecn);
    ~InclusionType();
    
    
public:
    
    std::string ITName;     ///< Type name
    int ITid;               ///< Type ID
    int ITN;                ///< In-plane symmetry
    double ITk;             ///< Type Kappa
    double ITkg;            ///< KappaG
    double ITk1;            ///< K_||
    double ITk2;            ///< K_norm
    double ITc0;            ///< Curvature
    double ITc1;            ///< Direction curvature 1
    double ITc2;            ///< Direction curvature 2
    double ITelambda;       ///< Line tension
    double ITekg;           ///< Geodesic rigidity
    double ITekn;           ///< Normal curvature line rigidity
    double ITecn;           ///< Spontaneous normal curvature
};

#endif
