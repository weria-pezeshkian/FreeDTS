#if !defined(AFX_ConstantExternalFieldOnVectorFields_H_334B21B8_INCLUDED_)
#define AFX_ConstantExternalFieldOnVectorFields_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "vertex.h"
#include "AbstractExternalFieldOnVectorFields.h"

class ConstantExternalFieldOnVectorFields : public AbstractExternalFieldOnVectorFields {

public:
    ConstantExternalFieldOnVectorFields(std::string data_stream);
    ~ConstantExternalFieldOnVectorFields();

    double GetCouplingEnergy(int layer, VectorField* p_vf, vertex *pvertex);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "ConstantFieldOnVectorFields";}
    inline static std::string GetDefaultReadName() {return "ConstantFieldOnVectorFields";}
    
private:
    std::vector<Vec3D> m_vFieldDirection;
    std::vector<double> m_vFieldStrength;

};


#endif
