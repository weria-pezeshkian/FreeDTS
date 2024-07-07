
#include "VAHGlobalMeshProperties.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "State.h"
VAHGlobalMeshProperties::VAHGlobalMeshProperties() : 
                                                                  m_TotalVolume(0.0),
                                                                  m_TotalArea(0.0),
                                                                  m_TotalCurvature(0.0) {

}
VAHGlobalMeshProperties::~VAHGlobalMeshProperties() {
    
}
void VAHGlobalMeshProperties::Initialize(State* pState){
    m_pState = pState;
    return;
}
void VAHGlobalMeshProperties::CalculateAVertexRingContributionToGlobalVariables(vertex *p_vertex, double &vol, double &area, double& curvature){
    
    vol = 0.0;
    area = 0.0;
    curvature = 0.0;
    double C = p_vertex->GetP1Curvature() + p_vertex->GetP2Curvature();
    curvature = C * (p_vertex->GetArea());
    
    const std::vector<triangle *>& ring_triangles = p_vertex->GetVTraingleList();
    for (std::vector<triangle *>::const_iterator it = ring_triangles.begin() ; it != ring_triangles.end(); ++it){
        vol+=CalculateSingleTriangleVolume(*it);
        area += (*it)->GetArea();
    }
    
    const std::vector<vertex *>& ring_vertex = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::const_iterator it = ring_vertex.begin() ; it != ring_vertex.end(); ++it){
         C = (*it)->GetP1Curvature() + (*it)->GetP2Curvature();
        double A = (*it)->GetArea();
        curvature += C * A;
    }

    
    return;
}
double VAHGlobalMeshProperties::CalculateSingleTriangleVolume(triangle *pTriangle){
    
    if(m_pState->GetMesh()->GetHasCrossedPBC()){
        *(m_pState->GetTimeSeriesLog()) << "---> the system has crossed the PBC while volume is being calculated.";
        *(m_pState->GetTimeSeriesLog()) << " SOLUTION: Restart the simulation and center the system. Also, activate the command for centering the box.";

         exit(0);
    }
    
    double T_area = pTriangle->GetArea();
    Vec3D Normal_v = pTriangle->GetNormalVector();
    Vec3D Pos = pTriangle->GetV1()->GetPos();

    // Compute triangle volume
    return T_area * (Pos.dot(Normal_v, Pos)) / 3.0;
}

void VAHGlobalMeshProperties::CalculateALinkTrianglesContributionToGlobalVariables(links *p_link, double &vol, double &area, double& curvature) {
    vol = 0.0;
    area = 0.0;
        
    vol += CalculateSingleTriangleVolume(p_link->GetTriangle());
    area += p_link->GetTriangle()->GetArea();
    
    vol += CalculateSingleTriangleVolume(p_link->GetMirrorLink()->GetTriangle());
    area += p_link->GetMirrorLink()->GetTriangle()->GetArea();
    
    vertex *v1 = p_link->GetV1();
    vertex *v2 = p_link->GetV2();
    vertex *v3 = p_link->GetV3();
    vertex *v4 = p_link->GetMirrorLink()->GetV3();

    double C1 = v1->GetP1Curvature() + v1->GetP2Curvature();
    double C2 = v2->GetP1Curvature() + v2->GetP2Curvature();
    double C3 = v3->GetP1Curvature() + v3->GetP2Curvature();
    double C4 = v4->GetP1Curvature() + v4->GetP2Curvature();
    double A1 = v1->GetArea();
    double A2 = v2->GetArea();
    double A3 = v3->GetArea();
    double A4 = v4->GetArea();
    curvature = C1*A1 + C2*A2 + C3*A3 + C4*A4;
    

    return;
}
void VAHGlobalMeshProperties::CalculateBoxRescalingContributionToGlobalVariables(double lx, double ly, double lz, double& vol, double& area, double& curvature){
    
    vol = 0;
    area = 0;
    curvature = 0;
    
    std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
        vol += CalculateSingleTriangleVolume(*it);
        area += (*it)->GetArea();
    }
    
    std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {
        double area = (*it)->GetArea();
        double curv = (*it)->GetP1Curvature() + (*it)->GetP2Curvature();
        curvature += curv*area;
    }
    
    return;
}
