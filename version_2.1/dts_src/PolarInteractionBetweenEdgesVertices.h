#ifndef AFX_PolarInteractionBetweenEdgesVertices_H_334B21B8_INCLUDED_
#define AFX_PolarInteractionBetweenEdgesVertices_H_334B21B8_INCLUDED_

#include "vertex.h"
#include "AbstractNonbondedInteractionBetweenVertices.h"

class PolarInteractionBetweenEdgesVertices : public AbstractNonbondedInteractionBetweenVertices {
public:
    PolarInteractionBetweenEdgesVertices(State* pstate, std::string input_data);
    ~PolarInteractionBetweenEdgesVertices();

    double GetVertexNonBondedEnergy(vertex* pvertex) const override;
    std::string CurrentState() const override;

    void Initialize() override;
    inline std::string GetDerivedDefaultReadName() const override { return "PolarInteractionBetweenEdgesVertices"; }
    inline static std::string GetDefaultReadName() { return "PolarInteractionBetweenEdgesVertices"; }

private:
    std::string m_Input_Data;
    State *m_pState;
    
    std::vector<vertex*> m_pEdgeV;
    double m_EP;
    double m_R0;
    double m_R0_2;
    double m_DAngle;
    Vec3D *m_pBox;
    
    
    double CalculateNonbondedInteractionBetweenTwoVertices(vertex *, vertex * ) const;

};

#endif
