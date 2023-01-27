#if !defined(AFX_SpringPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_)
#define AFX_SpringPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object is to couple the system to a harmonic potential between two groups, and also change the ...
 */
class SpringPotentialBetweenTwoGroups
{
public:
    SpringPotentialBetweenTwoGroups();
    SpringPotentialBetweenTwoGroups(bool state);
    SpringPotentialBetweenTwoGroups(bool state,  double K, double l0, int t0, std::string group1,std::string group2,double nx,double ny,double nz);
    ~SpringPotentialBetweenTwoGroups();

       inline bool GetState()                           {return m_State;}
    inline double GetEnergy()                           {return m_Energy;}
    inline double GetForce()                           {return m_Force;}

public:
    
    void MakeGroups(std::vector<vertex *> Ver,std::string ndx);
    void CalculateEnergy(int step);
    void MovingVertex(vertex* v, Vec3D Dx);
    void RejectMovingVertex(vertex* v, Vec3D Dx);

    //=====
private:
  
    std::vector<int> ReadIndex(std::string ndx,std::string groupname);
    Vec3D COMVertexGroup(std::vector<vertex *>);
    
    
private:
    double m_K;
    bool m_State;
    std::string m_Group1Name;
    std::string m_Group2Name;

    std::vector<vertex *> m_pGroup1;
    std::vector<vertex *> m_pGroup2;
    
    std::vector<int> m_Group1NDX;
    std::vector<int> m_Group2NDX;
    Vec3D m_Direction;
    double m_L0;
    double m_LT;
    double m_T0;
    Vec3D m_Group1COG;
    Vec3D m_Group2COG;

       double m_Energy;
    double m_Force;


    





};


#endif
