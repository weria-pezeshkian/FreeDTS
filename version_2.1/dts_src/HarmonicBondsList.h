#if !defined(AFX_HarmonicBondsList_H_334B21B8_INCLUDED_)
#define AFX_HarmonicBondsList_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "vertex.h"
#include "bond.h"
#include "AbstractBondedPotentialBetweenVertices.h"
class State;
class HarmonicBondsList : public AbstractBondedPotentialBetweenVertices {

public:
    HarmonicBondsList(State *pState, std::string filename);
    ~HarmonicBondsList();
    void Initialize();
  //  double GetBondEnergyOfVertex(vertex *pvertex);
    double GetTotalEnergy();
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName()  {return "HarmonicBondsList";}
    inline static std::string GetDefaultReadName() {return "HarmonicBondsList";}
    
private:
    std::string m_FileName;
    std::vector<bond> m_AllBonds;

    State *m_pState;
};


#endif
