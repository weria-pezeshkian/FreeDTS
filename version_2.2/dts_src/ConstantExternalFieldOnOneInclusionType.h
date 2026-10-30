#if !defined(AFX_ConstantExternalFieldOnOneInclusionType_H_884B21B0_INCLUDED_)
#define AFX_ConstantExternalFieldOnOneInclusionType_H_884B21B0_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "vertex.h"
#include "AbstractExternalFieldOnInclusions.h"

class ConstantExternalFieldOnOneInclusionType : public AbstractExternalFieldOnInclusions {

public:
    ConstantExternalFieldOnOneInclusionType(std::string inc_type, double k, double x, double y, double z);
    ~ConstantExternalFieldOnOneInclusionType();

    double GetCouplingEnergy(vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "ConstantExternalFieldOnOneInclusionType";}
    inline static std::string GetDefaultReadName() {return "ConstantExternalFieldOnOneInclusionType";}
    
private:
    Vec3D m_FieldDirection;
    double m_FieldStrength;
    std::string m_IncType;
};


#endif
