#if !defined(AFX_AbstractForceonVertices_H)
#define AFX_AbstractForceonVertices_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for Force on Vertices by Inclusions.
========================================================
*/
class State;
class  AbstractForceonVertices {
public:
    AbstractForceonVertices(){
        
    }
    virtual ~ AbstractForceonVertices(){
        
    }
    virtual double Energy_of_Force(vertex *p, Vec3D dx) = 0;
    virtual void Initialize() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;
    inline static std::string GetBaseDefaultReadName() {return "ForceOnVertices";}
    
private:
    
};
//---- a class for no box change
class NoForceonVertices : public AbstractForceonVertices {
public:
    NoForceonVertices(){
        
    }
    ~NoForceonVertices(){
        
    }
    
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline static std::string GetDefaultReadName()  {return "No";}

    void Initialize(){
        return;
    }
    double Energy_of_Force(vertex *p, Vec3D dx){
        return 0;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};
#endif
