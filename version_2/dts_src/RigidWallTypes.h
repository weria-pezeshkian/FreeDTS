#if !defined(AFX_FreeDTSRigidWallType_H)
#define AFX_FreeDTSRigidWallType_H
#include <iostream>
#include "AbstractBoundary.h"
// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for rigid wall boundry
========================================================
*/
//---- a class for no box change
class State;
class TwoFlatParallelWall : public  AbstractBoundary {
public:
    TwoFlatParallelWall(State* pState, double thickness, char direction);
    ~TwoFlatParallelWall();
    inline std::string GetDerivedDefaultReadName() {return "TwoFlatParallelWall";}
    inline static std::string GetDefaultReadName() {return "TwoFlatParallelWall";}

    void Initialize();
    bool MoveHappensWithinTheBoundary(double x, double y, double z, vertex* v);
    std::string CurrentState();

private:
    State* m_pState;
    double m_HalfThickness;
    double m_MidPlane;
    int m_Element;  // 0 X, 1, Y, 2 Z
    std::string m_Direction;  // 0 X, 1, Y, 2 Z

    
};

class EllipsoidalShell : public  AbstractBoundary {
public:
    EllipsoidalShell(State* pState, double thickness, double r, double a, double b, double c);
    ~EllipsoidalShell();
    inline std::string GetDerivedDefaultReadName() {return "EllipsoidalShell";}
    inline static std::string GetDefaultReadName() {return "EllipsoidalShell";}

    void Initialize();
    bool MoveHappensWithinTheBoundary(double dx, double dy, double dz, vertex* v);
    std::string CurrentState();

private:
    State* m_pState;
    double m_HalfThickness;
    double m_R;
    double m_A;
    double m_B;
    double m_C;
    Vec3D m_COG;
    Vec3D m_ElipScale;

    
};

#endif
