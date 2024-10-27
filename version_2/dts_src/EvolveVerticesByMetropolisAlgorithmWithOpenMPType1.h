#ifndef EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_OPENMP_1_H
#define EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_OPENMP_1_H


#include "SimDef.h"
#include "MESH.h"
#include "Vec3D.h"
#include "AbstractVertexPositionIntegrator.h"
#include "AbstractSimulation.h"

class State;
class EvolveVerticesByMetropolisAlgorithmWithOpenMPType1 : public AbstractVertexPositionIntegrator {
public:
    EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState);
    EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(State *pState, double rate_surf, double rate_edge, double dr);
    ~EvolveVerticesByMetropolisAlgorithmWithOpenMPType1();
    void Initialize();
    bool EvolveOneStep(int step);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithmOpenMP";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithmOpenMP";}

private:
    bool EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp, double &en);
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
    int m_Total_ThreadsNo;
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;
    
private:
    bool do_Simulation(){
        std::cout<<" ---> error, 999o1o this should have been called \n";
        return false;
    }

};

//, public MESH, public AbstractSimulation
#endif
