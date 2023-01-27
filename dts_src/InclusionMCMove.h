#if !defined(AFX_InclusionMCMove_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_)
#define AFX_InclusionMCMove_H_8P4B21B8_C13C_5648_BF23_914095086215__INCLUDED_


#include "SimDef.h"
#include "vertex.h"
#include "Energy.h"
#include "Vec3D.h"
#include "Inclusion_Interaction_Map.h"
#include "RNG.h"
class State;
class InclusionMCMove
{
public:
    InclusionMCMove ();
    InclusionMCMove (State *pState);
	 ~InclusionMCMove();


   void MC_Move_AnInclusion(inclusion *pinc, RNG *);

         inline double   GetEnergyDifference()    const    {return m_EnergyDifference;}
         inline int   GetMoveValidity()        {return m_MoveValidity;}





public:


private:

Inclusion_Interaction_Map * m_pInt;
    State *m_pState;
    double m_Beta;
    inclusion * m_pInc;
    double * m_ptotenergy;
private:
void KawasakiMove(double, links *);
void RotationMove(double,double,double);

double m_EnergyDifference;
int m_MoveValidity;


private:
    std::vector <links> m_LIntEChange;         /// links that may change their interaction energy after a vertex move
    std::vector <links*> m_pLIntEChange;      /// pointer to links that may change their interaction energy after a vertex move


    





};


#endif
