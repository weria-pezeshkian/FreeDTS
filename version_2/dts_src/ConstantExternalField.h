#if !defined(AFX_ConstantExternalField_H_334B21B8_INCLUDED_)
#define AFX_ConstantExternalField_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "vertex.h"
#include "AbstractExternalFieldOnInclusions.h"

class ConstantExternalField : public AbstractExternalFieldOnInclusions {

public:
    ConstantExternalField(double k, double x, double y, double z);
    ~ConstantExternalField();

    double GetCouplingEnergy(vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "ConstantExternalField";}
    inline static std::string GetDefaultReadName() {return "ConstantExternalField";}
    
private:
    Vec3D m_FieldDirection;
    double m_FieldStrength;

};


#endif
