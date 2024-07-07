#if !defined(AFX_AbstractInclusionConversion_H)
#define AFX_AbstractInclusionConversion_H
#include <iostream>

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
class  AbstractInclusionConversion {
public:
    AbstractInclusionConversion(){
        m_NumberOfAttemptedMoves = 0;
        m_AcceptedMoves = 0;
    }
    virtual ~ AbstractInclusionConversion(){
        
    }
    virtual void Initialize(State *pstate) = 0;
    virtual bool Exchange(int step) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "InclusionConversion";}

    
    
    
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    
};
//---- a class for no box change
class NoInclusionConversion : public AbstractInclusionConversion {
public:
    NoInclusionConversion(){
        
    }
    ~NoInclusionConversion(){
        
    }
    inline std::string GetDerivedDefaultReadName() {return "No";}
    inline static std::string GetDefaultReadName() {return "No";}


    
    void Initialize(State *pstate){
        return;
    }
    
    bool Exchange(int step){
        
        return false;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};
#endif
