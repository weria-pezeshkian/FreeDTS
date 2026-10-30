#if !defined(AFX_FlatVertexSubstrate_H_334B21B8_INCLUDED_)
#define AFX_FlatVertexSubstrate_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 September 2024
*/
#include "vertex.h"
#include "AbstractVertexAdhesionToSubstrate.h"

class FlatVertexSubstrate : public AbstractVertexAdhesionToSubstrate {

public:
    FlatVertexSubstrate(std::string data);
    ~FlatVertexSubstrate();

    double GetCouplingEnergy(vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "FlatSubstrate";}
    inline static std::string GetDefaultReadName() {return "FlatSubstrate";}
    
private:
    double m_AdhesionStrength;
    double m_Z;
    double m_Depth;

};


#endif
