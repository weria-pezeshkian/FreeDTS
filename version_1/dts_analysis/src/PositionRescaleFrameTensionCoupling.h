#if !defined(AFX_PositionRescaleFrameTensionCoupling_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_PositionRescaleFrameTensionCoupling_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "GenerateCNTCells.h"
#include "CouplingtoFixedGlobalCurvature.h"
#include "SpringPotentialBetweenTwoGroups.h"
#include "CNTCell.h"
#include "triangle.h"
//#include "State.h"
class State;
class PositionRescaleFrameTensionCoupling
{
public:
    PositionRescaleFrameTensionCoupling();
    PositionRescaleFrameTensionCoupling(double sigmap,State *st);
	~PositionRescaleFrameTensionCoupling();

    inline bool GetCNTCondition()        {return m_UpdateCNT;}



public:

    bool MCMoveBoxChange(double dr, double * TotalEnergy, double temp, int step, GenerateCNTCells *pGenCNT, std::vector<vertex *>, std::vector<links *>,std::vector<triangle* > );

private:
    void CheckCNTSize();
    void CheckMinDistance();
    double DistanceSquardBetweenTwoVertices(vertex *,vertex *,Vec3D );
    void CheckLinkLength();
    void CheckFaceAngle();
    void PerformMove();
    void RejectMove();
    void AcceptMove();
    bool CheckFaceAngle(links * l);
private:
Inclusion_Interaction_Map * m_pInc;

std::vector<CNTCell *> m_pAllCNT;  
std::vector<vertex* > m_pAllVertex;
std::vector<links* > m_pAllLinks;
std::vector<triangle* > m_pAllTriangle;
GenerateCNTCells *m_pGenCNT;

    
std::vector<triangle> m_AllTriangle;   
std::vector<vertex > m_AllVertex;
std::vector<links > m_AllLinks;
std::vector<links > m_AllProjectedLinks;   




    Vec3D *m_pBox;
double m_SigmaP;
    double m_dr;
    double m_drx;
    double m_dry;
    double m_Lyx;
    bool m_UpdateCNT;
    bool m_Move;
    
    
    double m_oldLx;
    double m_oldLy;
    double m_newLx;
    double m_newLy;
    int m_step;

    double m_Lnox;
    double m_Lnoy;
    

     double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
    CouplingtoFixedGlobalCurvature *m_pCFGC;
    SpringPotentialBetweenTwoGroups *m_pSPBTG;
    
    double m_DetaR;
    double m_DeltaA;



};


#endif
