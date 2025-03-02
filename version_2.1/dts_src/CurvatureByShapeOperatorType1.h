#if !defined(AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractCurvature.h"
#include "State.h"

/*
 * @file Curvature.h
 * @brief Declaration of the CurvatureByShapeOperatorType1 class.
 *
 * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * @copyright Weria Pezeshkian
 *
 * @class CurvatureByShapeOperatorType1
 * @brief Class to obtain curvature of a single vertex using shape operator method, introduced in Phys. Rev. E 81, 041922 (2010)
 *
 * This class computes the curvature of a vertex based on the shape operator method.
 * It provides functions to initialize the curvature computation and update the curvature of surface and edge vertices.
 *
 * This class assumes that the area of each triangle and edge shape operatror  has been correctly calculated.
 * Edge Shape operators get update in the links class.
 */
class State;
class CurvatureByShapeOperatorType1 : public AbstractCurvature {
public:
    
    CurvatureByShapeOperatorType1(State* pstate);
	 ~CurvatureByShapeOperatorType1();

public:
    bool Initialize();
    bool UpdateSurfVertexCurvature(vertex *p);
    bool UpdateEdgeVertexCurvature(vertex *p);
    bool UpdateVertexCurvature(vertex *p);

    inline  std::string GetDerivedDefaultReadName()  {return "ShapeOperator_1";}
    inline static std::string GetDefaultReadName()   {return "ShapeOperator_1";}
    std::string CurrentState();

private:
    Tensor2 Householder(const Vec3D &N);
    Vec3D Calculate_Vertex_Normal(vertex *p, double &area);
    State *m_pState;
    Vec3D *m_pBox;
};


#endif
