#if !defined(AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_)
#define AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "Vec3D.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Triangle object: points to 3 vertices and has an area and normal vector
 */

class triangle
{
public:
    
	triangle(int id, vertex *v1, vertex *v2, vertex *v3);
	triangle(int id);
	 ~triangle();

        inline const int GetTriID()          const  {return m_ID;}
        inline vertex *GetV1()                  	{return m_V1;}
        inline vertex *GetV2()                  	{return m_V2;}
        inline vertex *GetV3()                  	{return m_V3;}
        inline double GetArea()                  	{return m_Area;}
    	inline Vec3D GetAreaVector()                {return m_AreaVector;}
        inline Vec3D GetNormalVector()              {return m_Normal;}
        inline bool GetRepresentation()             {return m_Representation;}

public:
void UpdateRepresentation(bool); 	/// this is for visulaization output and does not effect the simulation
void UpdateNormal_Area(Vec3D *Box);     // A function to recalcualte the area and normal when the one of its vertices moves
void UpdateVertex(vertex *v1,vertex *v2,vertex *v3); // If a link flips the triangle changes its vertices, this function do the job


private:
    vertex *m_V1;
    vertex *m_V2;
    vertex *m_V3;
  int m_ID;
  bool m_Representation;
  Vec3D m_Normal;
  Vec3D m_AreaVector;
  double m_Area;

};


#endif
