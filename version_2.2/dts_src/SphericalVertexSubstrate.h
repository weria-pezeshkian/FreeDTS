#if !defined(AFX_SphericalVertexSubstrate_H_334B21B8_INCLUDED_)
#define AFX_SphericalVertexSubstrate_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 September 2024
*/
#include "vertex.h"
#include "AbstractVertexAdhesionToSubstrate.h"

class SphericalVertexSubstrate : public AbstractVertexAdhesionToSubstrate {

public:
    SphericalVertexSubstrate(std::string data);
    ~SphericalVertexSubstrate();

    double GetCouplingEnergy(vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "SphericalSubstrate";}
    inline static std::string GetDefaultReadName() {return "SphericalSubstrate";}
    
    friend class NonequilibriumCommands; // Friendship declaration

private:
    Vec3D m_Center;
    double m_Radius;
    double m_AdhesionStrength;

};


#endif
