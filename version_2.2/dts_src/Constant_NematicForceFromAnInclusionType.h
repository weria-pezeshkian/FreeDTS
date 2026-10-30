#if !defined(AFX_Constant_NematicForceFromAnInclusionType_H_444B21B8_INCLUDED_)
#define AFX_Constant_NematicForceFromAnInclusionType_H_444B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractForceonVerticesfromInclusions.h"
class Constant_NematicForceFromAnInclusionType : public AbstractForceonVerticesfromInclusions{
public:
    Constant_NematicForceFromAnInclusionType(double f0, std::string inc_typename);
    ~Constant_NematicForceFromAnInclusionType();
    double Energy_of_Force(vertex *p, Vec3D dx);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName()  {return "ConstantNematicForceFromAnInclusionType";}
    inline static std::string GetDefaultReadName() {return "ConstantNematicForceFromAnInclusionType";}
    Vec3D Inclusion_Force(vertex *pv2);

private:
    Vec3D ActiveNematicForce_1(vertex *pv2, vertex *pv1);
    double m_F0;
    double m_ActiveEnergy;
    std::string m_IncType;



};


#endif
