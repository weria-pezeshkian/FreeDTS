#ifndef EVOLVE_VERTICES_BY_KineticMonteCarlo_H
#define EVOLVE_VERTICES_BY_KineticMonteCarlo_H


#include "SimDef.h"
#include "MESH.h"
#include "Vec3D.h"
#include "AbstractVertexPositionIntegrator.h"
#include "AbstractSimulation.h"

class State;
class EvolveVerticesByKineticMonteCarlo : public AbstractVertexPositionIntegrator {
public:
    EvolveVerticesByKineticMonteCarlo(State *pState);
    EvolveVerticesByKineticMonteCarlo(State *pState, double DV, double dt);
    ~EvolveVerticesByKineticMonteCarlo();
    void Initialize();
    bool EvolveOneStep(int step);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "KineticMonteCarlo";}
    inline static std::string GetDefaultReadName() {return "KineticMonteCarlo";}

private:
    bool EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp);
    bool VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2);
    bool CheckFacesAfterAVertexMove(vertex* p_vertex);
    std::vector<links*> GetEdgesWithInteractionChange(vertex* p_vertex);

    double  SystemEnergy();  // it is for bug finding only; slow function, this is for development time, should be deleted
    State *m_pState;
  //  int Convert2LocalVoxelIndex(int new_index, int old_index, int No_index);
private:
    std::vector<vertex*>&        m_pSurfV;
    std::vector<vertex*>&       m_pEdgeV;
    Vec3D *m_pBox;
    
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;
    
    double m_dt;
    double m_Dv;
    double m_V0;
    double m_Vel_Standard;

    Vec3D GaussianDistribution_Vec(Vec3D Vel, double std);

private:
    bool do_Simulation(){
        std::cout<<" ---> error, 999o1o this should have been called \n";
        return false;
    }

};

//, public MESH, public AbstractSimulation
#endif
