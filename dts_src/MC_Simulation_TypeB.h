#if !defined(AFX_MC_Simulation_TypeB_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_)
#define AFX_MC_Simulation_TypeB_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "State.h"

class MC_Simulation_TypeB
{
public:
    
	MC_Simulation_TypeB(State *state);
	 ~MC_Simulation_TypeB();

public:


private:
    double m_Beta;
    double m_minAngle, m_Lmin2, m_Lmax2;
    Vec3D * m_pBox;
    std::vector<vertex*>      m_pAllV;
    std::vector<triangle*>    m_pAllT;
    std::vector<links*>       m_pAllLinks;
    std::vector<links*>       m_pHalfLinks1;
    std::vector<links*>       m_pHalfLinks2;
    std::vector<inclusion*>   m_pInclusions;
    
private:
    void  CenterIntheBox();
    double  SystemEnergy(State *pState);
    
    
    // A set of functions to check if a mesh is good for mc with vertices move
    bool    CheckMesh(MESH *pMesh);
    double  CheckLengthBetweenTwoVertex(vertex* v1, vertex* v2, Vec3D *pBox);
    double  CheckFaceAngle(links * l, Vec3D *);
    Vec3D   CalculateNormal(vertex* v1 ,vertex* v2 ,vertex* v3,Vec3D *pBox);
    void ReadIndexFile(std::string indexfilename);



};


#endif
