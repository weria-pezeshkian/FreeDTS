#ifndef ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H
#define ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Vec3D.h"
#include "AbstractAlexanderMove.h"

class State;

class AlexanderMoveByMetropolisAlgorithm : public AbstractAlexanderMove {
public:
    AlexanderMoveByMetropolisAlgorithm(State *pState);
    AlexanderMoveByMetropolisAlgorithm(State *pState, double rate);
    ~AlexanderMoveByMetropolisAlgorithm();

    bool Initialize();
    bool EvolveOneStep(int step);
    inline  std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithm";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithm";}
    
private:
    bool FlipOneEdge(int step, links *pedge, double temp);
    bool EdgeCanBeFliped(links *pedge);
    std::vector<links*> GetEdgesWithInteractionChange(links *p_edge);
    bool CheckFacesAfterFlip(links* pedge);
    double SystemEnergy(); // For bug finding only; slow function (should be deleted in production code)
    std::string CurrentState();

    
private:
    State *m_pState;

    std::vector<links*>&   m_pSurfL;
    Vec3D *m_pBox;    
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;
    
    
};

#endif // ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H
