#if !defined(AFX_AbstractTotalAreaCoupling_H)
#define AFX_AbstractTotalAreaCoupling_H
#include <iostream>
#include "VAHGlobalMeshProperties.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for global Curvature and how energy should change
========================================================
*/
class  AbstractTotalAreaCoupling  {
public:
    AbstractTotalAreaCoupling(VAHGlobalMeshProperties *pVHA) : m_pVAH(pVHA),  m_TotalArea(pVHA->m_TotalArea),  m_CalculatedGlobalVariable(pVHA->m_CalculatedGlobalVariable) {
        
    }
    virtual ~ AbstractTotalAreaCoupling(){
        
    }
    virtual  void Initialize(State *pstate) = 0;
    virtual  void UpdateTotalArea(double oldarea, double newarea) = 0;
    virtual  double CalculateEnergyChange(double oldarea,  double newarea) = 0;
    virtual double GetCouplingEnergy() = 0;
    virtual double CalculateAreaofALinkTriangles(links *p_link) = 0;
    virtual double CalculateAreaOfAVertexRing(vertex * pVeretx) = 0 ;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "TotalAreaCoupling";}
protected:
    VAHGlobalMeshProperties *m_pVAH;
    double &m_TotalArea;
    bool &m_CalculatedGlobalVariable;
};
//---- default value for no coupling to area change. 
class NoTotalAreaCoupling : public  AbstractTotalAreaCoupling {
public:
    NoTotalAreaCoupling(VAHGlobalMeshProperties *VHA) : AbstractTotalAreaCoupling(VHA) {
        
    }
    ~NoTotalAreaCoupling(){
        
    }
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline std::string GetDefaultReadName()  {return "No";}
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
    
    
    
    void Initialize(State *pstate){
        return;
    }
    double CalculateEnergyChange(double oldarea,  double newarea){
        return 0;
    }
    void UpdateTotalArea(double oldarea, double newarea){
        return;
    }
    double GetCouplingEnergy(){
        return 0;
    }
    double CalculateAreaofALinkTriangles(links *p_link){
        return 0;
    }
    double CalculateAreaOfAVertexRing(vertex * pVeretx){
        return 0;
    }
};

#endif
