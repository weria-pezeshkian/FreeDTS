
#include <limits>
#include "links.h"
#include "vertex.h"

links::links(int id, vertex *v1, vertex *v2, triangle *t1) {
    m_IntEnergy = 0;
    m_OldIntEnergy = 0;
    m_T1 = t1;
    m_V1 = v1;
    m_V2 = v2;
    m_ID = id;
    m_mirorflag = false;
    m_Show = true;
    m_EdgeSize = 0;
    m_LinkType = 0;
    m_Be = 0;
    m_He = 0;
    m_Number_of_VectorField_Layers = 0;
}
links::links(int id) {
    m_IntEnergy = 0;
    m_OldIntEnergy = 0;
    m_ID = id;
    m_mirorflag = false;
    m_Show = true;
    m_EdgeSize = 0;
    m_LinkType = 0;
    m_Be = 0;
    m_He = 0;
    m_Number_of_VectorField_Layers = 0;
}

links::~links() {
    
}
void links::UpdateTriangle(triangle *v){
    m_T1 = v;
    return;
}
void links::UpdateV(vertex *v1,vertex *v2,vertex *v3){
    m_V3=v3;
    m_V2=v2;
    m_V1=v1;
    return;
}
void links::UpdateV3(vertex *v3){
    m_V3=v3;
    return;
}
void links::InitializeVFIntEnergy(int no_vf){
    //-- in each run, this function should run only once 
    for (int i = 0; i<no_vf; i++){
        m_VFieldIntEnergy.push_back(0);
    }
    m_Number_of_VectorField_Layers = no_vf;
    return;
}

double links::GetVFIntEnergy() {
    /*
     * @brief Calculates the total vector field interaction energy.
     *
     * This function iterates through the vector of interaction energies and
     * accumulates the total energy for this specific connection.
     *
     * @return The total vector field interaction energy.
     */
    
    if(m_Number_of_VectorField_Layers == 0){
        return 0;
    }
    
    double en = 0.0; // Initialize total energy to zero

    // Iterate through the vector using iterators and accumulate the total energy
    for (std::vector<double>::iterator it = m_VFieldIntEnergy.begin(); it != m_VFieldIntEnergy.end(); ++it) {
        en += *it;
    }

    return en; // Return the total interaction energy
}
double links::GetVFIntEnergy(int layer) {

    if( layer >= m_Number_of_VectorField_Layers){
        std::cout<<"---> iligal action \n";
        return 0;
    }
    return m_VFieldIntEnergy[layer];
}
void links::UpdateMirrorLink(links* v){
    m_mirorlink=v;
    return;
}
bool links::Copy_InteractionEnergy(){
    
    m_OldIntEnergy = m_IntEnergy;
    return true;
}
bool links::Copy_VFInteractionEnergy(){
    
    if(m_Number_of_VectorField_Layers == 0){
        return true;
    }
    m_OldVFieldIntEnergy = m_VFieldIntEnergy;
    return true;
}
bool links::Reverse_VFInteractionEnergy(){
    
    if(m_Number_of_VectorField_Layers == 0){
        return true;
    }
    m_VFieldIntEnergy = m_OldVFieldIntEnergy;
    if(m_mirorflag){
        m_mirorlink->UpdateVFIntEnergy(m_OldVFieldIntEnergy);
    }
    return true;
}
bool links::Reverse_InteractionEnergy(){
    
    m_IntEnergy = m_OldIntEnergy;
    if(m_mirorflag){
        m_mirorlink->UpdateIntEnergy(m_IntEnergy);
    }
    return true;
}
void links::UpdateEdgeVector(Vec3D e_vector, double size){
    m_EdgeVector = e_vector;
    m_EdgeSize = size;
    
    return;
}
void links::ConstantMesh_Copy(){
   
   // m_OldLinkSide = m_LinkSide;
    m_OldNormal = m_Normal;
    m_OldBe = m_Be;
    m_OldHe = m_He;
    m_OldIntEnergy = m_IntEnergy;
    m_OldEdgeVector = m_EdgeVector;
    m_OldEdgeSize = m_EdgeSize;
        
    return;
}
void links::ReverseConstantMesh_Copy(){
    
    m_Be = m_OldBe;
    m_He = m_OldHe;
    m_IntEnergy = m_OldIntEnergy;
    m_Normal = m_OldNormal;
    m_EdgeVector = m_OldEdgeVector;
    m_EdgeSize = m_OldEdgeSize;
    
    //---- update some for mirror
    if(m_mirorflag){
        m_mirorlink->PutShapeOperator(m_Be,m_He);
        m_mirorlink->UpdateIntEnergy(m_IntEnergy);
        m_mirorlink->PutNormal(m_Normal);
        m_mirorlink->UpdateEdgeVector(m_EdgeVector, m_EdgeSize);
    }
    
    return;
}
bool links::SetCopy(){           // Copies the key ellements into the old type
    m_OldT1 = m_T1;     //
    m_OldV1 = m_V1;
    m_OldV2 = m_V2;
    m_OldV3 = m_V3;
    m_Oldneighborlink1 = m_neighborlink1;
    m_Oldneighborlink2  = m_neighborlink2;
    m_Oldmirorflag = m_mirorflag;
    m_OldNormal = m_Normal;
    m_OldBe = m_Be;
    m_OldHe = m_He;
    m_OldIntEnergy = m_IntEnergy;
    m_OldEdgeVector = m_EdgeVector;
    m_OldEdgeSize = m_EdgeSize;
    m_OldLinkType = m_LinkType;
    
    if(m_mirorflag){
        m_Oldmirorlink = m_mirorlink;
    }


    return true;
}
bool links::Reverse2PreviousCopy(){           // reverse to the last copy and part of the mirror variables
    
    m_T1 = m_OldT1;
    m_V1 = m_OldV1;
    m_V2 = m_OldV2;
    m_V3 = m_OldV3;
    m_neighborlink1 = m_Oldneighborlink1;
    m_neighborlink2 = m_Oldneighborlink2;
    m_mirorflag = m_Oldmirorflag;
    m_Normal = m_OldNormal;
    m_Be = m_OldBe;
    m_He = m_OldHe;
    m_IntEnergy = m_OldIntEnergy;
    m_EdgeVector = m_OldEdgeVector;
    m_EdgeSize = m_OldEdgeSize;
    m_LinkType = m_OldLinkType;
    
    //---- update some for mirror
    if(m_mirorflag){
        m_mirorlink = m_Oldmirorlink;
        m_mirorlink->PutShapeOperator(m_Be,m_He);
        m_mirorlink->UpdateIntEnergy(m_IntEnergy);
        m_mirorlink->PutNormal(m_Normal);
        m_mirorlink->UpdateEdgeVector(m_EdgeVector, m_EdgeSize);
    }

    return true;
}
void links::UpdateNeighborLink1(links* v){
    m_neighborlink1 = v;
}
void links::UpdateNeighborLink2(links* v){
    m_neighborlink2 = v;
}
void links::UpdateVisualize(bool v){
    m_Show = v;
}
void links::UpdateMirrorFlag(bool v){
    m_mirorflag = v;
    return;
}
void links::UpdateNormal()
{
   if(m_mirorflag){
       
       Vec3D v2 = (m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D v1 = m_T1->GetNormalVector();
       m_Normal = v1 + v2;
       double norm = m_Normal.norm();
       
       if(norm == 0){
           std::cout<<"error 2022----> one of the normals has zero size; normal link cannot be defined  \n";
           exit(0);
       }
       m_Normal.normalize();
       m_mirorlink->PutNormal(m_Normal);
    }
    else {
        // this is an edge link
        std::cout<<" developer error: link type and id "<<m_LinkType<<"  "<<m_ID<<" \n";
        std::cout<<"error ----> normal vector for edge links has not been defined   \n";
        exit(0);
    }
    return;
}
void links::PutNormal(Vec3D n) {
    
    m_Normal = n;
    return;
}
void links::PutShapeOperator(Vec3D Be,double He)
{
    m_Be = Be;
    m_He = He;
    
}
void links::UpdateIntEnergy(double en)
{
    m_IntEnergy = en;
}
bool links::UpdateVFIntEnergy(int vf_id, double en){
    
    if(vf_id >= m_VFieldIntEnergy.size()){
        std::cout<<"---> error, unexpected 293944 \n";
        return false;
    }
    
    m_VFieldIntEnergy[vf_id] = en;
    
    return true;
}
void links::UpdateVFIntEnergy(std::vector<double> VF_EN){

    m_VFieldIntEnergy = VF_EN;    
    return;
}
void links::UpdateShapeOperator(Vec3D *pBox)
{
    UpdateEdgeVector(pBox);
   if(m_mirorflag) {
       
       UpdateNormal();
       Vec3D Re = m_EdgeVector;
       Re=Re*(1.0/m_EdgeSize);
       Vec3D Be=m_Normal*Re;

    //====== Finding the size of the Be vector to make it normaized; this should not be needed
    // just have it so in case.
       double size=Be.norm();
       if(size!=0)
       {
           size=1.0/size;
       }
       else
       {
           std::cout<<" error 7634---> this should not happen \n";
           exit(0);
       }
       Be=Be*size;
       Vec3D Nf1=(m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D Nf2=m_T1->GetNormalVector();

       
//=== this is different from the orginal paper; it is faster
//==========
       double sign=Re.dot(Nf1*Nf2,Re);
       double tangle=Re.dot(Nf1,Nf2);
       double He=0;

       if(tangle < 1)
       {
           if(sign > 0) {
               He=-m_EdgeSize * sqrt(0.5*(1.0-tangle));
           }
           else if(sign < 0) {
            He=m_EdgeSize * sqrt(0.5*(1.0-tangle));  //He=2*cos(m_Dihedral/2.0)*renorm;
           }
           else
           {
            He = 0;
           }
       }
       else if(tangle >= 1 && tangle < 1.01)
       {
           // in case some numerical probelm happens; 1.01 is too large however,
           He = 0;
   
       }
       else if(tangle > 1.01)
       {
           std::cout<<"error--->: somthing wrong with this link \n";
           exit(0);
       }
	
       m_Be = Be;
       m_He = He;
       m_mirorlink->PutShapeOperator(m_Be , m_He);


  }
  else {

  }
}
void links::UpdateEdgeVector(Vec3D *pBox)
{
        double x1=m_V1->GetVXPos();
        double y1=m_V1->GetVYPos();
        double z1=m_V1->GetVZPos();
        double x2=m_V2->GetVXPos();
        double y2=m_V2->GetVYPos();
        double z2=m_V2->GetVZPos();
        
        double dx1=x2-x1;
        if(fabs(dx1)>(*pBox)(0)/2.0)
        {
            if(dx1<0)
            dx1=(*pBox)(0)+dx1;
            else if(dx1>0)
            dx1=dx1-(*pBox)(0);
        }
        double dy1=y2-y1;
        if(fabs(dy1)>(*pBox)(1)/2.0)
        {
            if(dy1<0)
            dy1=(*pBox)(1)+dy1;
            else if(dy1>0)
            dy1=dy1-(*pBox)(1);
        }
        double dz1=z2-z1;
        if(fabs(dz1)>(*pBox)(2)/2.0)
        {
            if(dz1<0)
            dz1=(*pBox)(2)+dz1;
            else if(dz1>0)
            dz1=dz1-(*pBox)(2);
        }
      
        Vec3D Re(dx1,dy1,dz1);
        m_EdgeVector = Re;
        m_EdgeSize=(m_EdgeVector.norm());
    
        if(m_mirorflag) {
            
            m_mirorlink->PutEdgeVector(m_EdgeVector*(-1),m_EdgeSize);
        }

}
void links::PutEdgeVector(Vec3D v, double l)
{
    m_EdgeSize   = l;
    m_EdgeVector = v;
}
bool links::Reverse_Flip(links *p_edge){
    /** Sine 2024
     * @brief Reverse the action of Flip(links *pedge) function.
     *
     *
     * @param p_edge A pointer to the edge to be flipped.
     * @return true if the reverse was successfully flipped, false otherwise.
     */
    // Check if the edge has a mirror link; if not, it is an edge link and cannot be flipped
    if (!p_edge->GetMirrorFlag()) {
        std::cerr << "Error (developer): An edge link without a mirror never been fliped to be reversed.\n";
        return false;
    }
    
    // Get the mirror edge and neighboring links
    links *p_medge = p_edge->GetMirrorLink(); // Mirror link
    links *l2 = p_edge->GetNeighborLink1();
    links *l3 = p_edge->GetNeighborLink2();
    links *l4 = p_medge->GetNeighborLink1();
    links *l1 = p_medge->GetNeighborLink2();

    // Get the triangles associated with the edge and its mirror
    triangle *t1 = p_edge->GetTriangle();
    triangle *t2 = p_medge->GetTriangle();

    // Get the vertices involved in the flip
    vertex *v4 = p_edge->GetV1();
    vertex *v3 = p_edge->GetV2();
    vertex *v1 = p_edge->GetV3();
    vertex *v2 = p_medge->GetV3();
    
    // Update the vertices of the edge and its mirror
    p_edge->UpdateV(v1, v2, v3);
    p_medge->UpdateV(v2, v1, v4);
    
    // Update the triangles with the new vertices
    t1->UpdateVertex(v1, v2, v3);
    t2->UpdateVertex(v2, v1, v4);

    // Update neighbor relationships for the vertices
    v1->AddtoNeighbourVertex(v2);
    v2->AddtoNeighbourVertex(v1);
    v3->RemoveFromNeighbourVertex(v4);
    v4->RemoveFromNeighbourVertex(v3);

    // Update the link lists for the vertices
    v1->AddtoLinkList(p_edge);
    v2->AddtoLinkList(p_medge);
    v3->RemoveFromLinkList(p_medge);
    v4->RemoveFromLinkList(p_edge);

    // Update the triangle lists for the vertices
    v3->RemoveFromTraingleList(t2);
    v4->RemoveFromTraingleList(t1);
    v1->AddtoTraingleList(t2);
    v2->AddtoTraingleList(t1);
    
    // Update neighbor links for the edge and its mirror
    p_edge->UpdateNeighborLink1(l1);
    p_edge->UpdateNeighborLink2(l2);
    p_medge->UpdateNeighborLink1(l3);
    p_medge->UpdateNeighborLink2(l4);

    // Update the third vertex of the neighboring links
    l1->UpdateV3(v1);
    l2->UpdateV3(v2);
    l3->UpdateV3(v2);
    l4->UpdateV3(v1);

    // Update neighbor links for the neighboring links
    l1->UpdateNeighborLink1(l2);
    l1->UpdateNeighborLink2(p_edge);
    l2->UpdateNeighborLink1(p_edge);
    l2->UpdateNeighborLink2(l1);
    l3->UpdateNeighborLink1(l4);
    l3->UpdateNeighborLink2(p_medge);
    l4->UpdateNeighborLink1(p_medge);
    l4->UpdateNeighborLink2(l3);

    // Update the triangles associated with the neighboring links
    l1->UpdateTriangle(t1);
    l3->UpdateTriangle(t2);

    return true;
}
bool links::Flip(links *p_edge) {
    /** Sine 2024
     * @brief Flips the specified edge in the mesh.
     *
     * This function flips the given edge by updating the associated vertices,
     * links, and triangles. If the edge does not have a mirror link (indicating
     * it is an edge link), the function will return false as flipping is not possible.
     *
     * @param p_edge A pointer to the edge to be flipped.
     * @return true if the edge was successfully flipped, false otherwise.
     */
    
    // Check if the edge has a mirror link; if not, it is an edge link and cannot be flipped
    if (!p_edge->GetMirrorFlag()) {
        std::cerr << "Error (developer): An edge link without a mirror cannot be flipped.\n";
        return false;
    }

    // Get the mirror edge and neighboring links
    links *p_medge = p_edge->GetMirrorLink(); // Mirror link
    links *l1 = p_edge->GetNeighborLink1();
    links *l2 = p_edge->GetNeighborLink2();
    links *l3 = p_medge->GetNeighborLink1();
    links *l4 = p_medge->GetNeighborLink2();

    // Get the triangles associated with the edge and its mirror
    triangle *t1 = p_edge->GetTriangle();
    triangle *t2 = p_medge->GetTriangle();

    // Get the vertices involved in the flip
    vertex *v1 = p_edge->GetV1();
    vertex *v2 = p_edge->GetV2();
    vertex *v3 = p_edge->GetV3();
    vertex *v4 = p_medge->GetV3();

    // Update the vertices of the edge and its mirror
    p_edge->UpdateV(v4, v3, v1);
    p_medge->UpdateV(v3, v4, v2);

    // Update the triangles with the new vertices
    t1->UpdateVertex(v4, v3, v1);
    t2->UpdateVertex(v3, v4, v2);

    // Update neighbor relationships for the vertices
    v1->RemoveFromNeighbourVertex(v2);
    v2->RemoveFromNeighbourVertex(v1);
    v4->AddtoNeighbourVertex(v3);
    v3->AddtoNeighbourVertex(v4);

    // Update the link lists for the vertices
    v1->RemoveFromLinkList(p_edge);
    v2->RemoveFromLinkList(p_medge);
    v4->AddtoLinkList(p_edge);
    v3->AddtoLinkList(p_medge);

    // Update the triangle lists for the vertices
    v1->RemoveFromTraingleList(t2);
    v2->RemoveFromTraingleList(t1);
    v3->AddtoTraingleList(t2);
    v4->AddtoTraingleList(t1);

    // Update neighbor links for the edge and its mirror
    p_edge->UpdateNeighborLink1(l2);
    p_edge->UpdateNeighborLink2(l3);
    p_medge->UpdateNeighborLink1(l4);
    p_medge->UpdateNeighborLink2(l1);

    // Update the third vertex of the neighboring links
    l1->UpdateV3(v4);
    l2->UpdateV3(v4);
    l3->UpdateV3(v3);
    l4->UpdateV3(v3);

    // Update neighbor links for the neighboring links
    l1->UpdateNeighborLink1(p_medge);
    l1->UpdateNeighborLink2(l4);
    l2->UpdateNeighborLink1(l3);
    l2->UpdateNeighborLink2(p_edge);
    l3->UpdateNeighborLink1(p_edge);
    l3->UpdateNeighborLink2(l2);
    l4->UpdateNeighborLink1(l1);
    l4->UpdateNeighborLink2(p_medge);

    // Update the triangles associated with the neighboring links
    l1->UpdateTriangle(t2);
    l3->UpdateTriangle(t1);

    return true;
}
void links::Flip()
{
    
   if(m_mirorflag){
    
    triangle *T2 = m_mirorlink->GetTriangle();
    vertex  *V4 = m_mirorlink->GetV3();
    vertex  *v1 = m_V1;
    vertex  *v2 = m_V2;
    vertex  *v3 = m_V3;
    vertex  *v4 = V4;
    links *l1=m_neighborlink1;
    links *l2=m_neighborlink2;
    links *l3=m_mirorlink->GetNeighborLink1();
    links *l4=m_mirorlink->GetNeighborLink2();

       m_V1->RemoveFromNeighbourVertex(m_V2);
       m_V2->RemoveFromNeighbourVertex(m_V1);
       V4->AddtoNeighbourVertex(m_V3);
       m_V3->AddtoNeighbourVertex(V4);

       m_V1->RemoveFromLinkList(this);
       m_V2->RemoveFromLinkList(m_mirorlink);
       V4->AddtoLinkList(this);
       m_V3->AddtoLinkList(m_mirorlink);
       m_V1->RemoveFromTraingleList(T2);
       m_V2->RemoveFromTraingleList(m_T1);
       m_V3->AddtoTraingleList(T2);
       V4->AddtoTraingleList(m_T1);

       this->UpdateNeighborLink1(l2);
       this->UpdateNeighborLink2(l3);
       m_mirorlink->UpdateNeighborLink1(l4);
       m_mirorlink->UpdateNeighborLink2(l1);
       l1->UpdateV3(V4);
       l2->UpdateV3(V4);
       l3->UpdateV3(m_V3);
       l4->UpdateV3(m_V3);
       l1->UpdateNeighborLink1(m_mirorlink);
       l1->UpdateNeighborLink2(l4);
       l2->UpdateNeighborLink1(l3);
       l2->UpdateNeighborLink2(this);
    
    l3->UpdateNeighborLink1(this);
    l3->UpdateNeighborLink2(l2);
    
    l4->UpdateNeighborLink1(l1);
    l4->UpdateNeighborLink2(m_mirorlink);
    

    l1->UpdateTriangle(T2);
    l3->UpdateTriangle(m_T1);
    l4->UpdateTriangle(T2);
    l2->UpdateTriangle(m_T1);
 
  
    m_V1=v4;
    m_V2=v3;
    m_V3=v1;
    V4=v2;
    m_mirorlink->UpdateV(m_V2,m_V1,V4);



    
    int id2=T2->GetTriID();
    int id1=m_T1->GetTriID();
    triangle tm1(id1,m_V1,m_V2,m_V3);
    triangle tm2(id2,m_V2,m_V1,V4);
      *m_T1 = tm1;
   
      *(m_mirorlink->GetTriangle()) = tm2;
    
   
   }
    else
    {
        std::cout<<"error---> a link without a mirror, possibly an edge link, is asked to be flipped, such an action is not possible \n";
        exit(0);
    }
    
}
//calaculates the dot(n1,n2) of this edge trinagle with the face of the mirror. this should not be smaller
// the minangle value
bool links::CheckFaceAngleWithMirrorFace(double &minangle) {
    // Check if mirror link exists
    if (!m_mirorflag) {
        return true;
    }
    // get normal vectors of the faces
    Vec3D n1 = m_T1->GetNormalVector();
    Vec3D n2 = m_mirorlink->GetTriangle()->GetNormalVector();
    
    if(n1.dot(n1,n2)<minangle){
        return false;
    }
    
    // Calculate and return the dot product of the normal vectors
    return true;
}
//calaculates the dot(n1,n2) of this edge trinagle with the face of the next edge. this should not be smaller
// the minangle value
bool links::CheckFaceAngleWithNextEdgeFace(double &minangle){
    
    if (!m_neighborlink1->GetMirrorFlag()) {
        return true; // such triangle does not exist, so true
    }
    Vec3D n1 = m_T1->GetNormalVector();
    Vec3D n2 = m_neighborlink1->GetMirrorLink()->GetTriangle()->GetNormalVector();
    
    if(n1.dot(n1,n2)<minangle){
        return false;
    }
    
    // Calculate and return the dot product of the normal vectors
    return true;
}
double links::Cal_CotOppositeAngle() {
    /**
     * @brief Computes the cotangent of the interior angle opposite to an edge.
     *
     * This function calculates the cotangent of the angle between two vectors
     * without explicit normalization, making it more efficient.
     *
     * @return Cotangent of the angle.
     */

    // Compute vectors from V3 to V1 and V3 to V2
    Vec3D X31 = m_V1->GetPos() - m_V3->GetPos();
    Vec3D X32 = m_V2->GetPos() - m_V3->GetPos();

    // Compute squared magnitudes (avoiding unnecessary sqrt calls)
    double S1 = Vec3D::dot(X31, X31);  // ||X31||^2
    double S2 = Vec3D::dot(X32, X32);  // ||X32||^2
    double S21 = Vec3D::dot(X31, X32); // X31 · X32

    // Compute the denominator for the cotangent formula
    double denom = S1 * S2 - S21 * S21;

    // Avoid precision errors and division by zero
    if (denom <= 1e-12) {
        return std::numeric_limits<double>::max();  // Return a large value for near-zero sine cases
    }

    // Compute and return cotangent as cos(theta) / sin(theta)
    return S21 / sqrt(denom);
}




