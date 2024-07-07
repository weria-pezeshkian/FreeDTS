 #if !defined(AFX_inclusion_CPP_7F4A21C7_C13Q_8823_BF2E_124095086234__INCLUDED_)
#define AFX_inclusion_CPP_7F4A21C7_C13Q_8823_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include "inclusion.h"
#include "vertex.h"
#include "Nfunction.h"
/*InclusionType::InclusionType()
{
    ITName  = "inc";   // type name
    ITid = -1;
    ITN = 0;      // inplane symmetry
    ITk = 0;
    ITkg  = 0;  // kappaG
    ITk1 = 0;  // K_||
    ITk2 = 0;  // K_norm
    ITc0 = 0;
    ITc1 = 0;
    ITc2 = 0;
}
InclusionType::InclusionType(std::string name, int id, int N, double k, double kg, double k1, double k2, double c0, double c1, double c2)
{
    ITName  = name;   // type name
    ITid = id;
    ITN = N;      // inplane symmetry
    ITk = k;
    ITkg  = kg;  // kappaG
    ITk1 = k1;  // K_||
    ITk2 = k2;  // K_norm
    ITc0 = c0;
    ITc1 = c1;
    ITc2 = c2;
}*/
inclusion::inclusion(int id)
{
m_ID=id;
m_TypeID = 0;
}

inclusion::~inclusion()
{
    
}
void inclusion::UpdateInclusionType(InclusionType* t)
{
    m_InclusionType = t;
}
void inclusion::UpdateInclusionTypeID(int Typeid)
{
    m_TypeID = Typeid;
}
void inclusion::Updatevertex(vertex * v)
{
    m_pvertex = v;
}
void inclusion::UpdateLocalDirection(Vec3D v)
{
    m_LDirection = v;
}
#endif



