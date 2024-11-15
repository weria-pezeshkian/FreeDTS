#if !defined(AFX_Boundary_H)
#define AFX_Boundary_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for rigid wall boundry
========================================================
*/
class  AbstractBoundary {
public:
    AbstractBoundary(){
        
    }
    virtual ~ AbstractBoundary(){
        
    }
    virtual void Initialize() = 0;
    virtual bool MoveHappensWithinTheBoundary(double x, double y, double z, vertex* v) = 0;
    virtual std::string CurrentState() = 0;
    virtual std::string CurrentStateParameters() = 0;

    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "Boundary";}
    inline static std::string GetErrorMessage(std::string info) {
        return "---> error: unknown AlexanderMove type -- \n";
    }
    

};
//---- a class for no box change
class PBCBoundary : public  AbstractBoundary {
public:
    PBCBoundary(){
        
    }
    ~PBCBoundary(){
        
    }

    inline std::string GetDerivedDefaultReadName() {return "PBC";}
    void Initialize(){
     
        return;
    }
    bool MoveHappensWithinTheBoundary(double x, double y, double z, vertex* v){
        return true;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
    std::string CurrentStateParameters(){
        std::string state = " 0 ";
        return state;
    }
};

#endif
