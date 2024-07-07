
#ifndef AFX_VectorFieldsRotationByMetropolisAlgorithm_H_INCLUDED_
#define AFX_VectorFieldsRotationByMetropolisAlgorithm_H_INCLUDED_

#include "AbstractVectorFieldsRotationMove.h"
#include "SimDef.h"

class State;
class links;
/*
 * @file VectorFieldsRotationByMetropolisAlgorithm.h
 * @brief Declaration of the VectorFieldsRotationByMetropolisAlgorithm class.
 *
 * This file contains the declaration of the VectorFieldsRotationByMetropolisAlgorithm class,
 * which implements an inclusion pose update algorithm based on the Metropolis algorithm.
 *
 * @author Weria
 */


class VectorFieldsRotationByMetropolisAlgorithm : public AbstractVectorFieldsRotationMove {
public:

    VectorFieldsRotationByMetropolisAlgorithm(State *pState);
    VectorFieldsRotationByMetropolisAlgorithm(State *pState, double rate, double dr);
    VectorFieldsRotationByMetropolisAlgorithm(State *pState, double rate);
    ~VectorFieldsRotationByMetropolisAlgorithm();


    bool Initialize();
    bool EvolveOneStep(int step);
    inline std::string GetDerivedDefaultReadName() { return "MetropolisAlgorithm"; }
    inline static std::string GetDefaultReadName() { return "MetropolisAlgorithm"; }
    std::string CurrentState();



private:

    bool RotationMove(int layer, vertex *p_vertex, double dx,  double thermal);
    
private:
    std::vector<links*> GetEdgesWithInteractionChange(links* p_link);
    
private:
    State *m_pState; ///< Pointer to the simulation state.
    std::vector<vertex*>&        m_pActiveV;
    double &m_Beta;
    double &m_DBeta;
    int &m_No_VectorFields_Per_V;

};

#endif
