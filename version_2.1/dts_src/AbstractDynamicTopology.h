#if !defined(AFX_DYNAMICTOPOLOGY_H)
#define AFX_DYNAMICTOPOLOGY_H

#include "SimDef.h"
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class vertex;
class AbstractDynamicTopology  {
public:
    AbstractDynamicTopology() {
        m_Surface_Genus = 0;
    };
    virtual ~AbstractDynamicTopology(){
        
    }

    virtual bool MCMove(int step) = 0;
    inline int  GetSurfaceGenus() {return m_Surface_Genus;}
    virtual void Initialize() = 0;
    virtual std::string CurrentState() = 0;

    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "DynamicTopology";}
    inline static std::string GetErrorMessage(std::string s) {
        return "---> error: unknown topology change type -- \n";
    }
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    double m_NumberOfMovePerStep;   // how many updates should be made per step
    int m_Surface_Genus;

};
//---- a class for no box change
class ConstantTopology : public AbstractDynamicTopology {
public:
    ConstantTopology() {
        
    }
    ~ConstantTopology(){
        
    }

    inline std::string GetDerivedDefaultReadName() {return "No";}
    inline static std::string GetDefaultReadName() {return "No";}
    void Initialize(){
        return;
    }
    bool MCMove(int step){
        return false;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};

#endif
