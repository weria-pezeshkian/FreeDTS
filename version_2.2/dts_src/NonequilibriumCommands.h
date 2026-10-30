#if !defined(AFX_NonequilibriumCommands_H)
#define AFX_NonequilibriumCommands_H
#include <iostream>
#include <vector>       // Required for std::vector
#include <functional>   // Required for std::function
#include <string>       // Required for std::string// A class to change the simulation and system variable on run time. 
#include "Vec3D.h"

/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
NonequilibriumCommands
========================================================
*/
class State;
class NonequilibriumCommands {


public:
    NonequilibriumCommands(State* pState);
    ~NonequilibriumCommands();

    void Run(int step);
    bool LoadCommand(std::string strcommand);

    
    
private:
    void ExpandEllipsoidalCoreWall(int rate, double dr);
    void ChangeTemperatureWithConstantRate(int rate, double DT);
    void ThinningEllipsoidalShell(int rate, double dr);
    void IncrementHarmonicPotentialBetweenTwoGroups(int Rate, double Dr);
    void IncrementVolumeCouplingSecondOrder(int Rate, double Dr);
    void IncrementSphericalSubstrateCenter(int Rate, Vec3D Dirct);

private:
    State* m_pState;
    int m_ActiveSimStep;

    // Replace member function pointers with std::function<void()>
    std::vector<std::function<void()>> m_FunctionContainer;
    std::vector<int> m_FunctionExecutionRate;
    std::vector<std::vector<std::string> > m_FunctionArguments;

};

#endif
