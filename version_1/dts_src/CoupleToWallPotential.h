#if !defined(AFX_CoupleToWallPotential_H_7D4B21B8_C13C_5648_BF23_444095086239__INCLUDED_)
#define AFX_CoupleToWallPotential_H_7D4B21B8_C13C_5648_BF23_444095086239__INCLUDED_


#include "SimDef.h"
#include "vertex.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 A class to couple the system to walls, for type of all exist
 TwoFlatParallelWall
 Cuboid
 Ellipsoid
 EllipsoidalShell
 This class makes sure that the vertex moves are within the boundry of the wall. Only TwoFlatParallelWall works with fixed tension simulation (dynamical box simulations).
 */
class CoupleToWallPotential
{
public:
    
        CoupleToWallPotential();
        CoupleToWallPotential(bool state, std::string info);
        ~CoupleToWallPotential();
        inline bool       GetState()                   {return m_State;}

public:
    void Initialize(std::vector <vertex *> Apv); // To find initial guess of the boundry and force the system within the eq time to reach the targeted boundry
    bool CheckVertexMoveWithinWalls(int step, double x, double y, double z, vertex* v); // at every mc step, it makes sure that the move is within the boundries

private:
    bool m_State;
    Vec3D m_COG;
    std::string m_PotentialType;
    std::vector<std::string> m_Data;
    int m_EQTime;
    
    // All the target variables. We evatually want to make the wall this after some moves
    double m_H;         // thickness of any shell
    double m_Lx;        // side length of a wall box, x componenet
    double m_Ly;
    double m_Lz;
    double m_A;         // 
    double m_B;
    double m_C;
    int m_Step;
    // All the current variables. 
    double m_h;         // active thickness of any shell
    double m_lx;        // side length of a wall box, x componenet
    double m_ly;
    double m_lz;
    double m_r;         //
    std::vector <vertex *> m_pAllV;
    bool m_ReachTargetWall;   // a variable to check if the system is within the targeted wall not active wall
    
private:
    bool AllVerticesAreInsideTheBound();
    bool APointIsInsideTheBound(Vec3D X);
    double DistanceOfAPointFromBound(Vec3D X);
    void MoveTheWallsTowardTheTarget(int step);
};


#endif
