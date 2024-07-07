#if !defined(AFX_VertexMCMove_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_)
#define AFX_VertexMCMove_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_


#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Energy.h"
#include "CNTCell.h"
#include "Vec3D.h"
#include "Inclusion_Interaction_Map.h"
#include "CouplingtoFixedGlobalCurvature.h"
#include "SpringPotentialBetweenTwoGroups.h"
#include "Curvature.h"
class State;
class VertexMCMove
{
public:
    VertexMCMove();
    VertexMCMove(State *pState);
    ~VertexMCMove();
    void MC_MoveAVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp,Curvature* mpcurv);

         inline double   GetEnergyDifference()    const    {return m_EnergyDifference;}
         inline int   GetMoveValidity()        {return m_MoveValidity;}

public:
    void   EnergyDifference(Curvature* mpcurv);
    void   AccpetMove();
    void   Move();
    void   RejectMove();
    bool   CheckDistnace();
    int   ValidMove();
//  void AddtoVertexList(vertex * z);

private:
Inclusion_Interaction_Map * m_pInc;
bool CheckLengthBetweenTwoVertex( vertex*);
void   UppdateVertexCNTCell();
bool   CheckFaceAngle(links *);

private:
double m_dx;
double m_dy;
double m_dz;
double m_X;
double m_Y;
double m_Z;
    double m_oldX;
    double m_oldY;
    double m_oldZ;
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
vertex *m_pvertex;
vertex  m_vertex;
std::vector <vertex> m_AVer;
std::vector<vertex *> m_pAVer;
std::vector <links> m_RingLinks;
std::vector <links> m_nLinks;
std::vector <triangle> m_Triangle;
double m_oldEnergy;
Vec3D *m_pBox;
bool m_ChangeCNT;
CNTCell* m_NewCNT;
CNTCell* m_OldCNT;
double *m_pTotEnergy;
double m_Thermal;
    double m_Beta;
double m_EnergyDifference;
    State *m_pState;
int m_MoveValidity;
    bool m_face;
    double m_NewVOL;
    int m_step;
std::vector <links> m_LIntEChange;         /// links that may change their interaction energy after a vertex move
std::vector <links*> m_pLIntEChange;      /// pointer to links that may change their interaction energy after a vertex move
    CouplingtoFixedGlobalCurvature *m_pCFGC;
    SpringPotentialBetweenTwoGroups *m_pSPBTG;
    double m_DetaR;
    double m_DeltaA;
    double m_simplexarea;
    double m_simplexvolume;
//  int m_ID;



    





};


#endif
