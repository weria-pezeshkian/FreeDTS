#if !defined(AFX_Constant_NematicForce_H_334B21B8_INCLUDED_)
#define AFX_Constant_NematicForce_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractForceonVerticesfromInclusions.h"
class Constant_NematicForce : public AbstractForceonVerticesfromInclusions{
public:
    Constant_NematicForce(double f0);
    ~Constant_NematicForce();
    double Energy_of_Force(vertex *p, Vec3D dx);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName()  {return "Constant_NematicForce";}
    inline static std::string GetDefaultReadName() {return "Constant_NematicForce";}
    
private:
    Vec3D ActiveNematicForce_1(vertex *pv2, vertex *pv1);
    double m_F0;
    double m_ActiveEnergy;

};


#endif
