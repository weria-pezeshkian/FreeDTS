#if !defined(AFX_AbstractExternalFieldOnVectorFields_H)
#define AFX_AbstractExternalFieldOnVectorFields_H
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
class  AbstractExternalFieldOnVectorFields {
public:
    AbstractExternalFieldOnVectorFields(){
        
    }
    virtual ~ AbstractExternalFieldOnVectorFields(){
        
    }
    virtual double GetCouplingEnergy(int layer, VectorField* p_vf, vertex *pvertex) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "ExternalFieldOnVectorFields";}

    
private:
    
};
//---- a class for no box change
class NoExternalFieldOnVectorFields : public AbstractExternalFieldOnVectorFields {
public:
    NoExternalFieldOnVectorFields(){
        
    }
    ~NoExternalFieldOnVectorFields(){
        
    }
    inline std::string GetDerivedDefaultReadName()  {return "No";};
    inline static std::string GetDefaultReadName() {return "No";}

    double GetCouplingEnergy(int layer, VectorField* p_vf, vertex *pvertex){
        return 0;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};
#endif
