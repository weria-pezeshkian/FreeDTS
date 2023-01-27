#if !defined(AFX_Analysis_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_)
#define AFX_Analysis_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "State.h"

#include "Nfunction.h"
#include "Vec3D.h"
#include "RNG.h"
#include "GenerateCNTCells.h"
#include "Curvature.h"
#include "Energy.h"
#include "LinkFlipMC.h"
#include "WritevtuFiles.h"
#include "Restart.h"
#include "BTSFile.h"
#include "Traj_XXX.h"
#include "CoupleToWallPotential.h"

class Analysis
{
public:
    
	Analysis(State *state);
	 ~Analysis();

public:

    
private:
    std::vector<MeshBluePrint> m_TRJ;
    void RemovePBC(MESH&);
    void Center(MESH&);
    void UpdateEnergy(MESH&);
    void UpdateVolume(MESH&);
    void UpdateEnergy(std::vector<triangle*> T);
    void Translate(std::vector<vertex*> v, Vec3D L);
    State *m_pState;
    double m_totEn;
    double m_totVolume;
    double m_totArea;
    void SeparateInOut(MESH&);
    void SeparateInOutByVer(MESH&);
    void Inclusion_order(MESH &mesh,double &total_H, double &no_inc, double &total_order);
    double Geo_Theta(vertex *v1, vertex *v2);



};


#endif
