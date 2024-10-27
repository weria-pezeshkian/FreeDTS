#if !defined(AFX_AbstractVertexAdhesionToSubstrate_H)
#define AFX_AbstractVertexAdhesionToSubstrate_H
#include <iostream>
#include "vertex.h"

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for curvature calculations.
========================================================
*/
class State;
class  AbstractVertexAdhesionToSubstrate {
public:
    AbstractVertexAdhesionToSubstrate(){
        
    }
    virtual ~ AbstractVertexAdhesionToSubstrate(){
        
    }
    virtual double GetCouplingEnergy(vertex *pvertex) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "VertexAdhesionToSubstrate";}

    
private:
    
};
//---- a class for no box change
class NoVertexAdhesionCoupling : public AbstractVertexAdhesionToSubstrate {
public:
    NoVertexAdhesionCoupling(){
        
    }
    ~NoVertexAdhesionCoupling(){
        
    }
    inline std::string GetDerivedDefaultReadName()  {return "No";};
    double GetCouplingEnergy(vertex *pvertex){
        return 0;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};
#endif
