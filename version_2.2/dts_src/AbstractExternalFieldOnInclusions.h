#if !defined(AFX_AbstractExternalFieldOnInclusions_H)
#define AFX_AbstractExternalFieldOnInclusions_H
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
class  AbstractExternalFieldOnInclusions {
public:
    AbstractExternalFieldOnInclusions(){
        
    }
    virtual ~ AbstractExternalFieldOnInclusions(){
        
    }
    virtual double GetCouplingEnergy(vertex *pvertex) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "ExternalFieldOnInclusions";}
    inline static std::string GetRegistryError(std::string type) {
        std::string txt = "---> error: unknown method for External Field type: " + type;
        return txt;
    }
    inline static std::string GetErrorMessage(std::string s) {
        return "---> error: unknown External Field type -- \n";
    }

    
private:
    
};
//---- a class for no box change
class NoExternalFieldOnInclusions : public AbstractExternalFieldOnInclusions {
public:
    NoExternalFieldOnInclusions(){
        
    }
    ~NoExternalFieldOnInclusions(){
        
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
