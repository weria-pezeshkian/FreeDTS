#if !defined(AFX_FlatInclusionSubstrate_H_334B21B8_INCLUDED_)
#define AFX_FlatInclusionSubstrate_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 September 2024
*/
#include "vertex.h"
#include "AbstractVertexAdhesionToSubstrate.h"

class FlatInclusionSubstrate : public AbstractVertexAdhesionToSubstrate {

public:
    FlatInclusionSubstrate(std::string data);
    ~FlatInclusionSubstrate();

    double GetCouplingEnergy(vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "FlatInclusionSubstrate";}
    inline static std::string GetDefaultReadName() {return "FlatInclusionSubstrate";}
    
private:
    double m_AdhesionStrength;
    double m_Z;
    std::string m_inctype;

};


#endif
