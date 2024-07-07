#if !defined(AFX_AbstractVectorFieldsRotationMove_H)
#define AFX_AbstractVectorFieldsRotationMove_H
#include <iostream>
#include "SimDef.h"

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
class  AbstractVectorFieldsRotationMove {
public:
    AbstractVectorFieldsRotationMove(){
        m_NumberOfMovePerStep = 1;
        m_DR = 0.1;
    }
    virtual ~ AbstractVectorFieldsRotationMove(){
        
    }
    virtual bool Initialize() = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "VectorFieldsRotation";}
    


    
    void SetMoveRate(double rate){
        m_NumberOfMovePerStep = rate;

        return;
    }
    double GetAcceptanceRate(bool reset) {
        // Check if there have been attempted moves to avoid division by zero
        if (m_NumberOfAttemptedMoves == 0) {
            return 0.0; // Return 0 if no moves have been attempted
        }

        // Calculate acceptance rate
        double rate = static_cast<double>(m_AcceptedMoves) / m_NumberOfAttemptedMoves;

        // Reset counters if requested
        if (reset) {
            m_NumberOfAttemptedMoves = 0;
            m_AcceptedMoves = 0;
        }

        return rate;
    }
    
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    double m_NumberOfMovePerStep;   // how many updates should be made per step
    double m_DR;

};

#endif
