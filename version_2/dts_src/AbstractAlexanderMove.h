#if !defined(AFX_AbstractAlexanderMove_H)
#define AFX_AbstractAlexanderMove_H
#include <iostream>

// Define a base class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for AlexanderMove
========================================================
*/
class State;
class  AbstractAlexanderMove {
public:
    AbstractAlexanderMove(){
        m_NumberOfMovePerStep = 1;
    }
    virtual ~ AbstractAlexanderMove(){
        
    }
    virtual bool Initialize() = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;

    
    
    inline static std::string GetBaseDefaultReadName() {return "AlexanderMove";}
    inline static std::string GetErrorMessage(std::string s) {
        return "---> error: unknown AlexanderMove type -- \n";
    }
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

};

#endif
