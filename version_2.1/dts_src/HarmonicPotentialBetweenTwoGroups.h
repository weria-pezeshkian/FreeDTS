#if !defined(AFX_HarmonicPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_)
#define AFX_HarmonicPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_
#include "SimDef.h"
#include "MESH.h"
#include "AbstractApplyConstraintBetweenGroups.h"
#include "NonequilibriumCommands.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object is to couple the system to a harmonic potential between two groups, and also change the ...
 */
class State;
class HarmonicPotentialBetweenTwoGroups : public AbstractApplyConstraintBetweenGroups {
public:
    HarmonicPotentialBetweenTwoGroups(State* pState, double K, double l0, std::string group1,std::string group2,double nx,double ny,double nz);
    ~HarmonicPotentialBetweenTwoGroups();

public:
    
    bool Initialize();
    double CalculateEnergyChange(vertex* p_vertex, Vec3D Dx);
    double CalculateEnergyChange(double lx, double ly, double z);
    void AcceptMove();
    std::string Output_DataString();
    inline  std::string GetDerivedDefaultReadName() {return "HarmonicPotentialBetweenTwoGroups";}
    inline static std::string GetDefaultReadName() {return "HarmonicPotentialBetweenTwoGroups";}
    std::string CurrentState();
    
    
    friend class NonequilibriumCommands; // Friendship declaration

    //=====
private:  
    Vec3D COMVertexGroup(std::vector<vertex *>);
    bool DistanceCheck(); // Checks if the distance between the centers of geomtry (COG) of two groups exceeds half the box size in any specified direction.
    
    
private:
    double m_K;
    std::string m_Group1Name;
    std::string m_Group2Name;
    std::vector<vertex *> m_pGroup1;
    std::vector<vertex *> m_pGroup2;
    double m_G1Size;
    double m_G2Size;
    Vec3D m_Direction;
    double m_L0;
    Vec3D m_Group1COG;
    Vec3D m_Group2COG;
    Vec3D m_T_Group1COG;  // a temparpty copy of m_Group1COG before accepting the move
    Vec3D m_T_Group2COG;
    double m_DE;
    State *m_pState;
    Vec3D *m_pBox;




};


#endif
