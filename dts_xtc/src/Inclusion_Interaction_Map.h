#if !defined(AFX_Inclusion_Interaction_Map_H_9P4B21B8_C13C_5648_BF23_124095086237__INCLUDED_)
#define AFX_Inclusion_Interaction_Map_H_9P4B21B8_C13C_5648_BF23_124095086237__INCLUDED_


#include "SimDef.h"
#include "Vec3D.h"
#include "vertex.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An object to read pair interaction info of the inclusion type interactions
 This object should be optimised later might even be changed .....
 */
struct PairInt {
    int FunctionType;
    std::vector <double> Varibale;
    
} ;
class Inclusion_Interaction_Map
{
public:

	Inclusion_Interaction_Map(std::string inputfilename);
    Inclusion_Interaction_Map();
    ~Inclusion_Interaction_Map();
public:
    PairInt GetPairInt (int,int);

private:
std::vector<std::vector<PairInt> > m_MAT;
};


#endif
