#if !defined(AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An class to obtain curvature of a single vertex
 This class only give the correct answer if the area of each triangle has been calculated correclty.
 */
class Curvature
{
public:
    
	Curvature();
	 ~Curvature();

public:
    void CalculateCurvature(vertex *p);
private:
    vertex * m_pVertex;
private:
    Tensor2 Householder(Vec3D N);
};


#endif
