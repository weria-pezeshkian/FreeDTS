
#include <stdio.h>
#include "InclusionType.h"
#include "Nfunction.h"
InclusionType::InclusionType() {
    
    ITName  =  "--";        // Type name
    ITid  = -1;                  // Type ID
    ITN  = 0 ;                   // In-plane symmetry
    ITk  = 0 ;                // Kappa
    ITkg = 0;               // KappaG
    ITk1 = 0;               // K_||
    ITk2 = 0;               // K_norm
    ITc0 = 0;               // Curvature
    ITc1 = 0;               // Direction curvature 1
    ITc2 = 0;               // Direction curvature 2
    ITelambda = 0;          // Line tension
    ITekg = 0;              // Geodesic rigidity
    ITekn = 0;              // Normal curvature line rigidity
    ITecn = 0;              // Spontaneous normal curvature
}
InclusionType::~InclusionType(){
    
}
InclusionType::InclusionType(std::string name, int id, int n, double k, double kg, double k1, double k2, double c0, double c1, double c2, double lam, double ekg, double ekn, double ecn){
    ITName  =  name;        // Type name
    ITid  = id;                  // Type ID
    ITN  = n ;                   // In-plane symmetry
    ITk  = k ;                // Kappa
    ITkg = kg;               // KappaG
    ITk1 = k1;               // K_||
    ITk2 = k2;               // K_norm
    ITc0 = c0;               // Curvature
    ITc1 = c1;               // Direction curvature 1
    ITc2 = c2;               // Direction curvature 2
    ITelambda = lam;          // Line tension
    ITekg = ekg;              // Geodesic rigidity
    ITekn = ekn;              // Normal curvature line rigidity
    ITecn = ecn;              // Spontaneous normal curvature
}

