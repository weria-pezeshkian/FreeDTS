#if !defined(AFX_Constant_NematicForceByVectorFields_H_334B21B8_INCLUDED_)
#define AFX_Constant_NematicForceByVectorFields_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractForceonVerticesfromVectorFields.h"
class Constant_NematicForceByVectorFields : public AbstractForceonVerticesfromVectorFields{
public:
    Constant_NematicForceByVectorFields(std::string data_stream);
    ~Constant_NematicForceByVectorFields();
    double Energy_of_Force(vertex *p, Vec3D dx);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName()  {return "Constant_NematicForce";}
    inline static std::string GetDefaultReadName() {return "Constant_NematicForce";}
    Vec3D VectorFields_Force(vertex *pv);

private:
    Vec3D ActiveNematicForce_1(int layer, vertex *pv2, vertex *pv1);

    std::vector<double> m_F0;
    double m_ActiveEnergy;

};


#endif
