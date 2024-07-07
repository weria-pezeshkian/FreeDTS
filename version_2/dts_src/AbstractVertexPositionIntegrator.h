#if !defined(AFX_AbstractVertexPositionIntegrator_H)
#define AFX_AbstractVertexPositionIntegrator_H

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
class  AbstractVertexPositionIntegrator {
public:
    AbstractVertexPositionIntegrator(){
        m_FreezGroupName = "";
        m_DR = 0.05;
        m_NumberOfMovePerStep_Surf = 1;
        m_NumberOfMovePerStep_Edge = 1;
    }
    virtual ~ AbstractVertexPositionIntegrator(){
        
    }
    virtual void Initialize() = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;

    
    inline std::string GetFreezGroupName()              {return m_FreezGroupName;}
    inline double GetDR()                               {return m_DR;}
    
    inline static std::string GetBaseDefaultReadName() {return "VertexPositionIntegrator";}
    
    //---- update
    void UpdateFreezGroupName(std::string name_freezgroup){
        m_FreezGroupName = name_freezgroup;
        return;
    }
    void UpdateDR(double dr){
        m_DR = dr;
        return;
    }
    void SetMoveRate(double update_rateSur, double update_rateEdge ){
        m_NumberOfMovePerStep_Surf = update_rateSur;
        m_NumberOfMovePerStep_Edge = update_rateEdge;
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
    std::string m_FreezGroupName;
    double m_DR;                // size of the box change of a vertex
    double m_NumberOfMovePerStep_Surf;   // how many updates should be made per step
    double m_NumberOfMovePerStep_Edge;   // how many updates should be made per step

};

#endif
