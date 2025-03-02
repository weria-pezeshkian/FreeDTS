#ifndef ABSTRACT_GLOBAL_CURVATURE_H
#define ABSTRACT_GLOBAL_CURVATURE_H

#include <iostream>
#include "VAHGlobalMeshProperties.h"

/*
 * AbstractGlobalCurvature: A base class for global curvature and energy changes.
 * Developed in 2024 by Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */
class  AbstractGlobalCurvature  {
public:
    AbstractGlobalCurvature(VAHGlobalMeshProperties *pVHA) : m_pVAH(pVHA), m_TotalCurvature(pVHA->m_TotalCurvature), m_TotalArea(pVHA->m_TotalArea),  m_CalculatedGlobalVariable(pVHA->m_CalculatedGlobalVariable) {
        
    }
    virtual ~ AbstractGlobalCurvature(){
        
    }
    
    virtual  void Initialize(State* pState) = 0;
    virtual double GetCouplingEnergy() = 0;
    virtual  double CalculateEnergyChange(double delta_area, double delta_curvature) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "GlobalCurvatureCoupling";}
    inline static std::string GetErrorMessage(std::string s) {
        return "---> error: unknown global curvature coupling type -- \n";
    }


protected:
    VAHGlobalMeshProperties *m_pVAH;
    double &m_TotalArea;
    double &m_TotalCurvature; // Delta A = h * m_TotalCurvature = h * Sum [2H_v * A_v]
    bool &m_CalculatedGlobalVariable;
    

};
//---- a class for no box change
class NoGlobalCurvature : public  AbstractGlobalCurvature {
public:
    NoGlobalCurvature(VAHGlobalMeshProperties *VHA) : AbstractGlobalCurvature(VHA) {}
    ~NoGlobalCurvature(){ }
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline static std::string GetDefaultReadName() {return "No";}
    void Initialize(State *pstate){return;}
    double CalculateEnergyChange(double delta_area, double delta_curvature){return 0;}
    double GetCouplingEnergy() {return 0;}

    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};

#endif
