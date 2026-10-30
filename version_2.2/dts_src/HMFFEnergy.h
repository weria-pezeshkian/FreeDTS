/*
  Compute the energy of a whole system in a scalar field using HMFF.

  Copyright (c) 2024-2025 European Molecular Biology Laboratory

  Author: Valentin Maurer <valentin.maurer@embl-hamburg.de>
*/
#if !defined(AFX_Energy_H_334B21B8_C13C_2248_BF23_347859023452__INCLUDED_)
#define AFX_Energy_H_334B21B8_C13C_2248_BF23_347859023452__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"
#include "Energy.h"

// this should be changed and I need to make declearation class 
#include "VertexInScalarFieldPotential.cpp"
#include "MRCParser.h"

class State;
class HMFFEnergy : public Energy {
public:
  HMFFEnergy(State *pState, std::string);

  ~HMFFEnergy(){
      delete m_pCalculator;
    }

public:
  inline std::string GetDerivedDefaultReadName() override { return "FreeDTS1.0_HMFF"; }
  inline static std::string GetDefaultReadName() { return "FreeDTS1.0_HMFF"; }

  double SingleVertexEnergy(vertex *p) override;
  Vec3D CalculateGradient(vertex *p_vertex);

private:
  HarmonicPotentialCalculator *m_pCalculator = nullptr;
  double m_CoherenceLambda = 0.0;
};

#endif
