#if !defined(AFX_LinkFlipMC_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_)
#define AFX_LinkFlipMC_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_


#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Energy.h"
#include "Vec3D.h"
#include "Inclusion_Interaction_Map.h"
#include "SpringPotentialBetweenTwoGroups.h"
#include "CouplingtoFixedGlobalCurvature.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 A class to perform link flip by MC move, 
 */
class State;
class LinkFlipMC
{
public:
    
    LinkFlipMC();
	LinkFlipMC(State *pState);
	 ~LinkFlipMC();

                inline double   GetEnergyDifference()    const    {return m_EnergyDifference;}
         	    inline int   GetMoveValidity()        {return m_MoveValidity;}





public:
    void   MC_FlipALink(int step, links *plinks,  double temp);
    void   EnergyDifference();
    void   AccpetMove();
    void   PerformMove();
    void   RejectMove();
    bool   CheckFlipCondition();
    bool   CheckFaceAngle();

//  void AddtoVertexList(vertex * z);

private:
    double m_Beta;
    int m_step;
Inclusion_Interaction_Map * m_pInc;
links *m_pLinks;
double m_oldEnergy;
    State *m_pState;
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
Vec3D *m_pBox;
std::vector <vertex> m_AVer;
std::vector <triangle> m_AT;
std::vector <links> m_AL;
std::vector <links> m_NeighborLinks;
std::vector <links> m_LIntEChange;         /// links that may change their interaction energy after a vertex move
std::vector <links*> m_pLIntEChange;      /// pointer to links that may change their interaction energy after a vertex move
    vertex *m_V1;
    vertex *m_V2;
    vertex *m_V3;
    vertex *m_V4;
    
    triangle *m_T1;
    triangle *m_T2;
    
    
    links *m_Mirror;
    links *m_L1;
    links *m_L2;
    links *m_L3;
    links *m_L4;
double *m_pTotEnergy;
double m_Thermal;
double m_EnergyDifference;
int m_MoveValidity;
    bool m_face;
    
    double m_OldVertexVolume;
    double m_NewVertexVolume;
    double m_DV;
    CouplingtoFixedGlobalCurvature *m_pCFGC;
    SpringPotentialBetweenTwoGroups *m_pSPBTG;

    double m_DetaR;
    double m_DeltaA;

    double m_simplexarea;
    double m_simplexvolume;


};


#endif
