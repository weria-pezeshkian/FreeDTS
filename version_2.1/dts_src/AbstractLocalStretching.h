#if !defined(AFX_AbstractLocalStretching_H)
#define AFX_AbstractLocalStretching_H
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
class  AbstractLocalStretching {
public:
    AbstractLocalStretching(){
        
    }
    virtual ~ AbstractLocalStretching(){
        
    }
    virtual double Energy(vertex *p) = 0;
    virtual void Initialize() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;
    inline static std::string GetBaseDefaultReadName() {return "LocalStretching";}
    
private:
    
};
//---- a class for no box change
class NoLocalStretching : public AbstractLocalStretching {
public:
    NoLocalStretching(){
        
    }
    ~NoLocalStretching(){
        
    }
    
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline static std::string GetDefaultReadName()  {return "No";}

    void Initialize(){
        return;
    }
    double Energy(vertex *p){
        return 0;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};
#endif
