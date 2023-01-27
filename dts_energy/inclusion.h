#if !defined(AFX_inclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_inclusion_H_8P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_


#include "SimDef.h"
#include "CNTCell.h"
#include "Vec3D.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 inclusion object and type
 */
struct InclusionType {
    std::string ITName ;   // type name
    int ITid;       // Type ID
    int ITN;      // inplane symmetry
    double ITk;     // Type Kappa
    double ITkg;  // kappaG
    double ITk1;  // K_||
    double ITk2;  // K_norm
    double ITc0;  // curvature
    double ITc1;    // direction curvature 1
    double ITc2;    // direction curvature 2
};
class vertex;
class inclusion
{
public:
    
	inclusion(int id);
	 ~inclusion();

	    inline const int GetID()                               const  {return m_ID;}
        inline vertex* Getvertex()                                    {return m_pvertex;}
        inline Vec3D GetLDirection()                                  {return m_LDirection;}
        inline InclusionType* GetInclusionType()                       {return m_InclusionType;}
        inline int GetInclusionTypeID()                                {return m_TypeID;}

public:
    
  void UpdateInclusionType(InclusionType* );
  void UpdateInclusionTypeID(int Typeid);
  void Updatevertex(vertex * );
  void UpdateLocalDirection(Vec3D );

public:
    void ReadInclusionFromFile(std::ifstream *inputfile,std::vector <vertex *> pv);
    void WriteInclusionToFile(std::ofstream *output);
    void WriteInclusion();


private:
    int m_TypeID;
    int m_ID;
    InclusionType *m_InclusionType;
    Vec3D m_LDirection;      /// its direction in the local frame
    Vec3D m_GDirection;       /// its direction in the global frame
    vertex *m_pvertex;


};


#endif
