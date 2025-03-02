


#include "CurvatureByShapeOperatorType1.h"
#include "Tensor2.h"
CurvatureByShapeOperatorType1::CurvatureByShapeOperatorType1(State *pState) : m_pState(pState) {
    
}
CurvatureByShapeOperatorType1::~CurvatureByShapeOperatorType1(){

}
bool CurvatureByShapeOperatorType1::Initialize(){
 /*
  Curvature Calculation Initialization:

  This function initializes the curvature and area on all the vertices and edges. The sequence involves updating various geometric properties of the mesh components before computing curvature values.
  
     1) Update triangle area
     2) Update edge shape operator
     3) update vertex shape operator
    

    Note: it also initalize m_pState and m_pBox for later use.
    To compute the curvature of any vertex accurately, we must first compute the area and normal vectors of all vertices, along with the properties of all edges within a single ring (excluding the ring perimeter).
  */

        m_pBox = m_pState->GetMesh()->GetBox();
    
    // update all the triangles area and normal;
    const std::vector<triangle *>& pAllTriangles = m_pState->GetMesh()->GetActiveT();
    for (std::vector<triangle *>::const_iterator it = pAllTriangles.begin(); it != pAllTriangles.end(); ++it) {
        
        (*it)->UpdateNormal_Area(m_pBox);
    }
//---> update the normal of the edges and then update the shape operator
    const std::vector<links *>& pAllRight_edges = m_pState->GetMesh()->GetRightL();
    for (std::vector<links *>::const_iterator it = pAllRight_edges.begin(); it != pAllRight_edges.end(); ++it) {
        
        //(*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    //---> edge links
        const std::vector<links *>& pAlllink_edges = m_pState->GetMesh()->GetEdgeL();
        for (std::vector<links *>::const_iterator it = pAlllink_edges.begin(); it != pAlllink_edges.end(); ++it) {
            
                (*it)->UpdateEdgeVector(m_pBox);
        }
//---> now the surface vertices curvature, including vertex shape operator
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetSurfV();
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        
            UpdateSurfVertexCurvature(*it);
    }
//---> update edge vertices curvature
    const std::vector<vertex *>& pAlledge_vertex = m_pState->GetMesh()->GetEdgeV();
    for (std::vector<vertex *>::const_iterator it = pAlledge_vertex.begin(); it != pAlledge_vertex.end(); ++it) {
        UpdateEdgeVertexCurvature(*it);
    }
    
    return true;
}
bool CurvatureByShapeOperatorType1::UpdateVertexCurvature(vertex *pvertex){
 
    if(pvertex->GetVertexType() == 0){
        return UpdateSurfVertexCurvature(pvertex);
    }
    else if(pvertex->GetVertexType() == 1){
        return UpdateEdgeVertexCurvature(pvertex);
    }
    else{
        std::cout<<" vertex type "<<pvertex->GetVertexType()<<" is unrecognized \n";
        return false;
    }
    
    return true;
}
bool CurvatureByShapeOperatorType1::UpdateSurfVertexCurvature(vertex * pvertex){
    
    // Get vertex area from Calculate_Vertex_Normal
    double Area;
    
    // Update vertex normal and area
    Vec3D Normal = Calculate_Vertex_Normal(pvertex, Area);
    
    //--- update this values in the vertex
    pvertex->UpdateNormal_Area(Normal,Area);


    Tensor2  SV;
    Tensor2 IT('I');
    Tensor2 P=IT-IT.makeTen(Normal);
    Tensor2 Pt=P.Transpose();
    std::vector<links *> NLinks=pvertex->GetVLinkList();

    // Assuming necessary includes and context are already provided

    // Check if the edge vertex has one triangle
    for (auto it = NLinks.begin(); it != NLinks.end(); ++it) {
        links* currentLink = *it; // Store the pointer in a variable to avoid dereferencing multiple times
        if (currentLink->GetMirrorFlag()) {
            const Vec3D& ve = currentLink->GetNormal(); // Use reference to avoid copying
            double we = ve.dot(Normal, ve); // Calculate dot product
            
            const Vec3D& Be = currentLink->GetBe(); // Use reference to avoid copying
            double he = currentLink->GetHe(); // Get the height
            Vec3D Se = P * Be; // Matrix-vector multiplication
            
            double ff = Se.norm(); // Get the norm
            if (ff == 0) {
                std::cerr << "-----> Error: projection is zero error" << "\n";
                return false;
            } else {
                Se = Se*(1.0 / ff); // Normalize the vector
            }
            
            Tensor2 Q = P.makeTen(Se); // Create the tensor
            SV = SV + Q * (we * he); // Update SV
        } else {
            std::cerr << "---> error: vertex, with id " << pvertex->GetVID() << " is wrongly here \n";
            return false;
        }
    }



    //==== Compute Curvature and Local Frame

    Tensor2 Hous = Householder(Normal);
    Tensor2 LSV = Hous.Transpose() * (SV * Hous); // Local curvature matrix (2×2 minor matrix)

    // Compute eigenvalues of 2×2 symmetric matrix
    double a = LSV(0,0), b = LSV(0,1), c = LSV(1,1);
    double trace = a + c;
    double determinant = a * c - b * b;
    double delta = trace * trace - 4 * determinant;

    double k1, k2;
    if (delta >= 0.0) {
        // If delta is positive, compute eigenvalues normally
        double sqrt_delta = sqrt(delta);
        k1 = 0.5 * (trace + sqrt_delta); // Larger eigenvalue
        k2 = 0.5 * (trace - sqrt_delta); // Smaller eigenvalue
    } else {
        // Handle small numerical errors when delta is near zero
        k1 = k2 = 0.5 * trace;
        std::cerr << "---> WARNING: Failed to find curvature on vertex "
                  << pvertex->GetVID() << " due to negative delta (" << delta
                  << "). \n";
    }

    // Compute Eigenvectors for the 2×2 system
    Tensor2 EigenvMat('O');
    double v1 = b, v2 = k1 - a;
    double magnitude = sqrt(v1 * v1 + v2 * v2);

    // Ensure numerical stability in eigenvector computation
    if (magnitude < 1e-8) {
        v1 = 1.0;
        v2 = 0.0;
        magnitude = 1.0;
    }

    EigenvMat(0,0) =  v1 / magnitude;
    EigenvMat(1,0) =  v2 / magnitude;
    EigenvMat(0,1) = -EigenvMat(1,0);
    EigenvMat(1,1) =  EigenvMat(0,0);
    EigenvMat(2,2) =  1.0;

    // Compute coordinate transformation matrices
    Tensor2 TransferMatLG = Hous * EigenvMat;  // Local-to-Global
    Tensor2 TransferMatGL = TransferMatLG.Transpose();  // Global-to-Local

    pvertex->UpdateL2GTransferMatrix(TransferMatLG);
    pvertex->UpdateG2LTransferMatrix(TransferMatGL);

    // Normalize curvatures by area
    pvertex->UpdateP1Curvature(k1 / Area);
    pvertex->UpdateP2Curvature(k2 / Area);


    return true;
}
bool CurvatureByShapeOperatorType1::UpdateEdgeVertexCurvature(vertex * pvertex)
{
    // first we obtain the vertex area and normal. Area is not important here as the vertex is an edge vertex
    double Area=0.0;
    Vec3D Normal = Calculate_Vertex_Normal(pvertex, Area);
    pvertex->UpdateNormal_Area(Normal,Area);
    
    // the shape of the system                          //         v
                                                        //     l1 / \ l2
                                                 
    links* link1 = pvertex->m_pPrecedingEdgeLink;
    links* link2 = pvertex->m_pEdgeLink;

    double l1 = link1->m_EdgeSize;
    double l2 = link2->m_EdgeSize;
    double vlenght = 0.5*(l1+l2);
    pvertex->m_VLength = vlenght;
    Vec3D L1=link1->m_EdgeVector;
    Vec3D L2=link2->m_EdgeVector;
    L1 = L1*(1/l1);
    L2 = L2*(1/l2);
    Vec3D Norm = (L1-L2)*(1.0/vlenght);    // dT/ds = -Norm ; the size is curvature
    Vec3D Tv = Norm*Normal;  // T at each vertex
    Tv = Tv*(1/Tv.norm());
    Vec3D P = Normal*Tv;
    double cn = Norm.dot(Norm,Normal);
    double cg = Norm.dot(Norm,P);
    pvertex->m_Geodesic_Curvature = cg;
    pvertex->m_Normal_Curvature =cn;
    Tensor2 TransferMatGL(Tv,P,Normal);     // P1,P2,N is for other surfaces
    Tensor2 TransferMatLG=TransferMatGL.Transpose();
    pvertex->UpdateL2GTransferMatrix(TransferMatLG);
    pvertex->UpdateG2LTransferMatrix(TransferMatGL);
    
    return true;
    
}
Tensor2 CurvatureByShapeOperatorType1::Householder(const Vec3D& N){
    
    Vec3D Zk;
    Zk(2) = 1.0;
    Zk = Zk + N;
    Zk.normalize();
    Tensor2 I('I');
    Tensor2 W = Tensor2::makeTen(Zk);
    
    return (W*2-I);
}
Vec3D CurvatureByShapeOperatorType1::Calculate_Vertex_Normal(vertex *pvertex, double &area){
    // first we obtain the vertex area and normal.
    std::vector<triangle *> Ntr=pvertex->GetVTraingleList();
    
    area=0.0;
    Vec3D Normal;
    
    const std::vector<triangle *>& pNTriangles = pvertex->GetVTraingleList();
    for (auto* tri : pvertex->GetVTraingleList()) {
        const Vec3D& Nv = tri->GetAreaVector();
        Normal = Normal + Nv;
        area += tri->GetArea();
    }
    area = area/3.0;
    // Check for non-positive area
    if(area < 1e-8 ){
        
        *(m_pState->GetTimeSeriesLog())<<"---> error: vertex, with id "<<pvertex->GetVID() <<" has a negetive or zero area \n";
        return Normal;
    }
    double normalsize=Normal.norm();
    if(normalsize < 1e-8 ){
        
        *(m_pState->GetTimeSeriesLog())<<"---> error: vertex, with id "<<pvertex->GetVID() <<" has zero normal \n";
        return Normal;
    }
    Normal = Normal*(1.0/normalsize);
    
    return Normal;
}
std::string CurvatureByShapeOperatorType1::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}




