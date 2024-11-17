#if !defined(AFX_UserDefinedForceonVertices_H_334B21B8_INCLUDED_)
#define AFX_UserDefinedForceonVertices_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractForceonVertices.h"
class UserDefinedForceonVertices : public AbstractForceonVertices{
public:
    UserDefinedForceonVertices(State *pState, std::string inputs);
    ~UserDefinedForceonVertices();
    double Energy_of_Force(vertex *p, Vec3D dx);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName()  {return "User";}
    inline static std::string GetDefaultReadName() {return "User";}
    
    
    void Initialize();
private:
    Vec3D CalculateForce(vertex *pv);
    std::string m_Inputs;
    std::vector<vertex*> m_pallV;
    State *m_pState;

};


#endif