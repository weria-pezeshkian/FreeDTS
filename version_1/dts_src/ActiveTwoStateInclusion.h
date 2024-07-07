#if !defined(AFX_ActiveTwoStateInclusion_H_999B21B8_C13C_5648_BC23_444775086239__INCLUDED_)
#define AFX_ActiveTwoStateInclusion_H_999B21B8_C13C_5648_BC23_444775086239__INCLUDED_


//#include "SimDef.h"
#include "Inclusion_Interaction_Map.h"
#include "inclusion.h"
#include "RNG.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Active inclusion
 This command is not tested and may contain errors, do not used it .....
 */

class ActiveTwoStateInclusion
{
public:
    
    ActiveTwoStateInclusion();
    ActiveTwoStateInclusion(bool state, double ep1, double ep2, double per, double gama, std::string t1,std::string t2 );
	 ~ActiveTwoStateInclusion();



    inline bool  GetHealth()                    {return m_Health;}
    inline bool  GetState()                    {return m_State;}
    inline double  *GetActiveEnergy()           {return &m_ActiveEnergy;}
    inline int    GetDeltaN()                   {return m_DeltaN;}






public:
    void ActiveExchange(double *energy);
void Initialize(std::vector<inclusion*> pInc, std::vector<InclusionType*> inctype, Inclusion_Interaction_Map * pint, RNG* rng);


private:
    InclusionType *m_pIncType1;
    InclusionType *m_pIncType2;
    std::vector<InclusionType*> m_AllIncType;
    std::vector<inclusion*> m_pInclusions;
    bool m_Health;
    RNG *m_RNG;
    double m_Epsilon1;
    double m_Epsilon2;
    double m_Percentge;
    double m_Gama;
    std::string m_IncType1Name;
    std::string m_IncType2Name;
    std::vector<inclusion*> m_pInc;
    int m_N2;
    int m_N1;
    int m_DeltaN0;
    int m_N;
    int totDN;
    Inclusion_Interaction_Map *m_pInt;
    double m_ActiveEnergy;
    int m_DeltaN;
    bool m_State;



};


#endif
