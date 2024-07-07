#if !defined(AFX_AbstractApplyConstraintBetweenGroups_H)
#define AFX_AbstractApplyConstraintBetweenGroups_H

#include <iostream>
#include <string>

// Define an abstract class with virtual functions and some main functions
/*
=======================================================
 Developed 2024 by Weria Pezeshkian
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class is an abstract class for applying constraints between groups.
========================================================
*/

class State;  // Forward declaration of State class

class AbstractApplyConstraintBetweenGroups {
public:
    AbstractApplyConstraintBetweenGroups() {}
    virtual ~AbstractApplyConstraintBetweenGroups() {}

    inline double GetEnergy(){
        return m_Energy;
    }
    virtual bool Initialize() = 0;
    virtual double CalculateEnergyChange(vertex* p_vertex, Vec3D Dx) = 0;
    virtual double CalculateEnergyChange(double lx, double ly, double z) = 0;
    virtual void AcceptMove() = 0;

    virtual std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;

    static std::string GetBaseDefaultReadName() {
        return "ConstraintBetweenGroups";
    }

protected:
    double m_Energy;
};

// Class for no constraint
class NoConstraint : public AbstractApplyConstraintBetweenGroups {
public:
    NoConstraint() {}
    ~NoConstraint() {}

    std::string GetDerivedDefaultReadName() {return "No";}
    static std::string GetDefaultReadName() {return "No";}

    bool Initialize() {
        return true;
    }

    double CalculateEnergyChange(vertex* p_vertex, Vec3D Dx) {
        return 0.0;
    }
    double CalculateEnergyChange(double lx, double ly, double z){
        return 0.0;
    }
    void AcceptMove() {
        // No action needed for NoConstraint
    }
    void AcceptMove(double lx, double ly, double z){
        // No action needed for NoConstraint
    }

    std::string CurrentState() {
        return GetBaseDefaultReadName() + " = " + GetDerivedDefaultReadName();
    }
};

#endif // !defined(AFX_AbstractApplyConstraintBetweenGroups_H)
