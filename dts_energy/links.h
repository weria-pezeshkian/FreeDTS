#if !defined(AFX_links_H_6R4B21B8_K13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_links_H_6R4B21B8_K13C_5648_BF23_124095086234__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "Vec3D.h"
#include "Tensor2.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
link object (edge):
 */

class links
{
public:
    
	links(int id, vertex *v1, vertex *v2, triangle *t1);
	links(int id);
	 ~links();

// a set of inline function to access the private meembers
        inline const int GetID()                  const  {return m_ID;}
        inline vertex *GetV1()                       {return m_V1;}
        inline vertex *GetV2()                       {return m_V2;}
        inline vertex *GetV3()                       {return m_V3;}
        inline  bool GetVisualize ()                {return m_Show;}
        inline bool GetMirrorFlag ()           {return m_mirorflag;}
        inline triangle *GetTriangle()         {return m_T1;}
        inline links   *GetMirrorLink()                       {return m_mirorlink;}
        inline links   *GetNeighborLink1()                    {return m_neighborlink1;}
        inline links   *GetNeighborLink2()                    {return m_neighborlink2;}
        inline Vec3D GetNormal()                   {return m_Normal;}
    	inline Vec3D GetBe()                   	   {return m_Be;}
    	inline double GetHe()                      {return m_He;}
        inline int GetSimTimeStep()                      {return m_SimTimeStep;}
        inline int GetLinkSide()                      {return m_LinkSide;}
        inline double GetIntEnergy()                      {return m_IntEnergy;}

    


public:
void UpdateTriangle(triangle *v);
void UpdateV(vertex *v1,vertex *v2,vertex *v3);
void UpdateV3(vertex *v3);
void UpdateMirrorLink(links *v);
void UpdateNeighborLink1(links *v);
void UpdateNeighborLink2(links *v);
void UpdateVisualize(bool v);
void UpdateNormal();
void PutNormal(Vec3D);
void UpdateShapeOperator(Vec3D *);
void PutShapeOperator(Vec3D Be,double He);
void Flip();
void UpdateMirrorFlag(bool v);
void UpdateSimTimeStep(int v);
void UpdateIntEnergy(double en);

void ReadLinkFromFile(std::ifstream *inputfile,std::vector <vertex *> pv, std::vector <links *> pL, std::vector <triangle *> pT);
void WriteLinkToFile(std::ofstream *inputfile);
private:



private:

triangle *m_T1;
vertex  *m_V1;
vertex  *m_V2;
vertex  *m_V3;
links   *m_mirorlink;
links   *m_neighborlink1; /// triangle is 1-2-3;  the link is 1->2   this point to 2->3    
links   *m_neighborlink2; /// the link is 1->2   this is 3->1


    int m_ID;
    bool m_Show;
    bool m_mirorflag;
    int m_LinkSide;

    private:
    Vec3D m_Normal;
    
    Vec3D m_Be;
    double m_He;
    double m_IntEnergy;
    int m_SimTimeStep;



};


#endif
