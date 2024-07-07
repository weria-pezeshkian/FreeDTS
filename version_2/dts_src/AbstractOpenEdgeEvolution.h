#if !defined(AFX_AbstractOpenEdgeEvolution_H)
#define AFX_AbstractOpenEdgeEvolution_H
#include <iostream>
#include "RNG.h"

// Define a base class with a virtual function for different open edge treatment algorthems
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractOpenEdgeEvolution {
public:
    AbstractOpenEdgeEvolution(){
        m_NumberOfAttemptedMoves = 0;
        m_AcceptedMoves = 0;
        m_EdgeSize = 0;
    }
    virtual ~AbstractOpenEdgeEvolution(){
        
    }
    virtual bool Move(int step) = 0;
    virtual void Initialize() = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "OpenEdgeEvolution";}
    inline const int GetEdgeSize()            const {return m_EdgeSize;}  // note, this is number of edge links not the length

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

    
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    int m_EdgeSize;

};
//---- a class for no edge change
class NoEvolution : public AbstractOpenEdgeEvolution {
public:
    NoEvolution(){
        
    }
    ~NoEvolution(){
        
    }

    void Initialize(){
        return;
    }
    bool Move(int step){
        return false;
    }
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline static std::string GetDefaultReadName()  {return "No";}
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};

#endif
