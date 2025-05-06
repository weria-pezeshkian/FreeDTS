#ifndef AFX_InteractionBetweenInclusionsIn3D_H_334B21B8_INCLUDED_
#define AFX_InteractionBetweenInclusionsIn3D_H_334B21B8_INCLUDED_

#include "vertex.h"
#include "AbstractNonbondedInteractionBetweenVertices.h"

class InteractionBetweenInclusionsIn3D : public AbstractNonbondedInteractionBetweenVertices {
public:
    InteractionBetweenInclusionsIn3D(State* pstate, std::string input_data);
    ~InteractionBetweenInclusionsIn3D();

    double GetVertexNonBondedEnergy(vertex* pvertex) const override;
    std::string CurrentState() const override;

    void Initialize() override;
    inline std::string GetDerivedDefaultReadName() const override { return "InteractionBetweenInclusionsIn3D"; }
    inline static std::string GetDefaultReadName() { return "InteractionBetweenInclusionsIn3D"; }

private:
    
    
    std::vector<vertex*> FindCloseVertices(vertex* pV) const;

    std::string m_Input_Data;
    State *m_pState;
    
    double m_EP;

    
    
    double CalculateNonbondedInteractionBetweenTwoVertices(vertex *, vertex * ) const;

};

#endif
