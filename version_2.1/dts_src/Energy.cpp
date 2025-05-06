

#include <time.h>
#include "Energy.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Energy of a single vertex
 Energy of a link (when connected vertices has inclusions)
 Energy of the whole system
 */
Energy::Energy(State* pState)  : m_Box(pState->GetMesh()->Link2ReferenceBox()) {
    
    m_pState = pState;
}
Energy::~Energy() {
    
}
double Energy::SingleVertexEnergy(vertex *p_vertex) {
    
    double Energy=0.0;
//---> note, if the vertex has no inclusion, this will return zero by default
    Energy += m_pState->GetExternalFieldOnInclusions()->GetCouplingEnergy(p_vertex);
    Energy += m_pState->GetVertexAdhesionToSubstrate()->GetCouplingEnergy(p_vertex);

    Energy += m_pState->GetAbstractLocalStretching()->Energy(p_vertex);

    if(p_vertex->m_VertexType == 0) {
    
        Energy += SurfVertexBendingAndStretchingEnergy(p_vertex);
    }
    else if( p_vertex->m_VertexType == 1) {
        
        Energy += EdgeVertexBendingAndStretchingEnergy(p_vertex);
    }
    else{
        std::cout<<"---> error (developer) this is unexpected e6p5748 \n";
    }

    p_vertex->UpdateEnergy(Energy);

    return Energy;
}
double Energy::SurfVertexBendingAndStretchingEnergy(vertex * p_vertex){
    
    double en = 0;

    double c1 = p_vertex->GetP1Curvature();
    double c2 = p_vertex->GetP2Curvature();

    double mean_times2 = c1 + c2;
    double gussian = c1 * c2;
    double area = p_vertex->GetArea();
    
    // energy for area coupling K*(area-a0)^2
    if(m_Ka != 0){
        en += m_Ka*(area-m_Area0)*(area-m_Area0);
    }
    
    if(!p_vertex->VertexOwnInclusion()) {
        //---> if the vertex is just a membrane; kappa/2*(2H-c0)^2-kg*K
        en += (m_kappa*(mean_times2-m_SCurvature0)*(mean_times2-m_SCurvature0) - m_kappa_G * gussian) * area;
    }
    else{
        
        // e_v = kappa/2*(2H-c0)^2-kgK+k1/2(cp-cp0)^2+k2/2(cn-cn0)^2
        inclusion *p_inc = p_vertex->GetInclusion();
        double k0 = p_inc->m_IncType->ITk;
        double kg = p_inc->m_IncType->ITkg;
        double k1 = p_inc->m_IncType->ITk1;
        double k2 = p_inc->m_IncType->ITk2;
        double c0 = p_inc->m_IncType->ITc0;
        double cp10 = p_inc->m_IncType->ITc1;
        double cn20 = p_inc->m_IncType->ITc2;
        //--- kappa/2*(2H-c0)^2-kgK+
        double ev = k0 * (mean_times2-c0) * (mean_times2-c0) - kg * gussian;
        en += ev * area;
        
//--> if the k1 and k2 are zero, we do not need to calculate the rest
        if(k1 != 0 || k2 != 0){  // k1/2(cp-cp0)^2+k2/2(cn-cn0)^2
        
            Vec3D local_direction = p_inc->GetLDirection();
            double Cos = local_direction(0);
            double Sin = local_direction(1);
            double Cp = c1 * Cos * Cos + c2 * Sin * Sin;
            double Cn = c2 * Cos * Cos + c1 * Sin * Sin;
            double Delta_Cp = Cp - cp10;
            double Delta_Cn = Cn - cn20;
            en += (k1 * Delta_Cp * Delta_Cp + k2 * Delta_Cn * Delta_Cn) * area;
        }
    }
    
    return en;
}
double Energy::EdgeVertexBendingAndStretchingEnergy(vertex *p_vertex)
{
    
    double en = 0;

    double geo_c = p_vertex->GetGeodesicCurvature();
    double norm_c = p_vertex->GetNormalCurvature();
    double length = p_vertex->GetLength();
    
    // streaching energy for length coupling K*(l-l0)^2
    if(m_Kl != 0){
        en += m_Kl*(length-m_l0)*(length-m_l0);
    }
    
    if(!p_vertex->VertexOwnInclusion()) {
        // e = lammda + k * k_g*k_g + k2 * k_n*k_n
        double en_b = m_Lambda + m_Kappa_Geo * geo_c * geo_c + m_Kappa_Norm * norm_c * norm_c;
        
        en += en_b * length;
    }
    else {
        
        inclusion *p_inc = p_vertex->GetInclusion();
        double lambda = p_inc->m_IncType->ITelambda;
        double kg     = p_inc->m_IncType->ITekg;
        double kn     = p_inc->m_IncType->ITekn;
        double cn0    = p_inc->m_IncType->ITecn;
        
        en += (lambda + kg * geo_c * geo_c) * length;
        
        if(kn != 0) {
            
            Vec3D LD = p_inc->GetLDirection();
            double Cos = LD(0);
            kn = kn*Cos*Cos;
            double Delta_norm_c = (norm_c-cn0);
            
            en += (kn*Delta_norm_c*Delta_norm_c)*length;
        }
        else{  // I do not understand what does this means. 
            en += kn * norm_c * norm_c * length;
        }
    }


    return en;
}

double Energy::CalculateAllLocalEnergy()
{
    double E = 0.0;

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    const std::vector<links *>& pRight_L = m_pState->GetMesh()->GetRightL();
    const std::vector<links *>& pEdge_L = m_pState->GetMesh()->GetEdgeL();

    int number_of_vectorfields = m_pState->GetMesh()->GetNoVFPerVertex();
    
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
        
        E += SingleVertexEnergy(*it);
        //---- this is for the vector fields; in VertexVectorFields the vectorfield energy will be set. such a function does not set the binding energy in the vectorfield class
        E += CalculateVectorFieldMembraneBindingEnergy(*it);
    }
     for (std::vector<links *>::const_iterator it = pRight_L.begin() ; it != pRight_L.end(); ++it) {
         
        E += TwoInclusionsInteractionEnergy(*it);
         //-- interaction energies of vector fields
         for (int i = 0; i<number_of_vectorfields; i++){
             E += TwoVectorFieldInteractionEnergy(i, *it);
         }
    }
    for (std::vector<links *>::const_iterator it = pEdge_L.begin() ; it != pEdge_L.end(); ++it) {
        E += TwoInclusionsInteractionEnergy(*it);
        //-- interaction energies of vector fields
        for (int i = 0; i<number_of_vectorfields; i++){
            E += TwoVectorFieldInteractionEnergy(i, *it);
        }

    }
    return E;
}
double Energy::TwoInclusionsInteractionEnergy(links * p_edge) {

    vertex * p_v1 = p_edge->GetV1();
    vertex * p_v2 = p_edge->GetV2();

//--> If any of the vertices does not have an inclusion, the interaction is zero
    if(p_v1->VertexOwnInclusion()==false || p_v2->VertexOwnInclusion()==false){
       
        if(p_edge->GetMirrorFlag()==true) {
            p_edge->UpdateIntEnergy(0);
            (p_edge->GetMirrorLink())->UpdateIntEnergy(0);
        }
        else {
            p_edge->UpdateIntEnergy(0);
        }
        return 0;
    }

//--> now use the type if to get the interaction parameters of this pars
    
    int id1=(p_v1->GetInclusion())->m_IncType->ITid;
    int id2=(p_v2->GetInclusion())->m_IncType->ITid;
    PairInt pair_ab = m_pInt->GetPairInt(id1,id2);
    std::vector <double> ff = pair_ab.Varibale;
    //---> get the type of function that they interact, it is an intger
    int FunctionType  = pair_ab.FunctionType;

        double e_int = 0;
    switch (FunctionType) {
        case 0:{
            e_int = 0;
             break;
        }
         case 1: {
             double theta = (ff[2] != 0) ? Geo_Theta(p_v1, p_v2) : 0.0;
             e_int = InteractionFunction(ff[0], ff[1], ff[2], theta);

             break;
         }
        case 4: {
            double theta = (ff[2] != 0) ? Geo_Theta(p_v1, p_v2) : 0.0;
            e_int = InteractionFunctionFull(ff[0], ff[1], ff[2], theta, p_edge);

            break;
        }
        case 5: {
            e_int = InteractionFive(ff[0], ff[1], ff[2],ff[3],  p_edge);

            break;
        }
        case 3:{
              e_int = Filament_int(ff[0], ff[1], ff[2],p_v1, p_v2);
             break;
        }
        case 2:{
             m_Angle3D = 0;
             m_Angle2D = 0;
            e_int = F2(p_v1, p_v2, ff);
             break;
        }
        case 10:{
             m_Angle3D = 0;
             m_Angle2D = 0;
            e_int = F10(p_v1, p_v2, ff);
             break;
        }
        case 11:{
             m_Angle3D = 0;
             m_Angle2D = 0;
            e_int = F11(p_v1, p_v2, ff);
             break;
        }
        default:{
             std::cerr << "---> Error: Unrecognized function type ID --> " << FunctionType << std::endl;
             exit(0);
        }
     }
    
    if(p_edge->GetMirrorFlag()) {
        p_edge->UpdateIntEnergy(e_int/2.0);
        (p_edge->GetMirrorLink())->UpdateIntEnergy(e_int/2.0);
    }
    else {
        // divided by 2 here is also fine as everywhere will be muliplied by two again
        p_edge->UpdateIntEnergy(e_int/2.0);
    }
    return e_int;
}
double Energy::InteractionFunction(double N, double A, double B, double theta) {
    
        double e = 0;
        e = cos( double(N) * theta );
        e = -A + B * e;

    return e;
}
double Energy::InteractionFunctionFull(double N, double A, double B, double theta, links* pl) {
//---- to have the metric in
    double scale =  pl->Cal_CotOppositeAngle();
    if(pl->GetMirrorFlag()){
        scale += pl->GetMirrorLink()->Cal_CotOppositeAngle();
        scale  = scale/2;
    }
    
    double e =  -A + B * cos( double(N) * theta );

    return scale*e;
}
double Energy::InteractionFive(double N, double A, double B, double C,  links* p_link){
   
    vertex *v1 = p_link->GetV1();
    vertex *v2 = p_link->GetV2();

    Vec3D X1 = v1->GetPos();
    Vec3D X2 = v2->GetPos();
    Vec3D geodesic_dir=(X2-X1);

    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>  m_Box(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i) = geodesic_dir(i) + m_Box(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i) = geodesic_dir(i) - m_Box(i);
        }
    }
    double edge_size = geodesic_dir.norm();
    
    Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
    y1(2)=0;
    y1.normalize();
    Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
    y2(2)=0;
    y2.normalize();
    Vec3D n(0,0,1);
    
    Vec3D d1 = (v1->GetInclusion())->GetLDirection();
    Vec3D d2 = (v2->GetInclusion())->GetLDirection();
    double cos1 = y1.dot(y1,d1);
    double sin1 = n.dot(n*y1,d1);
    double cos2 = y1.dot(y2,d2);
    double sin2 = n.dot(n*y2,d2);

    // Compute theta
    double S_an = sin1 * sin2 + cos1 * cos2;
    double theta = acos(S_an);

    
    double e = cos( double(N) * theta );

    double cos2T_1 = cos1*cos1 + cos2*cos2 - 1; //2cos(theta1)^2-1+2cos(theta2)^2 -1, with below
    cos2T_1 = 2 * cos2T_1;
    
    e = -A + B * e + C * cos2T_1/edge_size;

return e;
    
}
double Energy::Filament_int(double A, double B, double C, vertex* p_v1, vertex* p_v2) {
    
        double theta = Geo_Theta(p_v1, p_v2);
        Vec3D d1 = (p_v1->GetL2GTransferMatrix())*(p_v1->GetInclusion()->GetLDirection());
        Vec3D d2 = (p_v2->GetL2GTransferMatrix())*(p_v2->GetInclusion()->GetLDirection());
        double ang_3d = Vec3D::dot(d1,d2);
        double e = -A + B * cos( 2 * theta ) + C*(1-ang_3d*ang_3d);

    return e;
}
double Energy::F2(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = -e0+e1*cosN(phi-phi0)+e2*(l/l0-1)cos(beta-beta0)
    /// 

    double E = 0;
    m_Angle2D = Geo_Theta(v1,v2);
    Vec3D N1 = v1->GetNormalVector();
    Vec3D N2 = v2->GetNormalVector();
    double beta = acos(N1.dot(N1,N2));
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double theta0 = var.at(3);
    double e2 = var.at(4);
    double beta0 = var.at(5);
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";

    theta0 = theta0/180.0*3.14;
    beta0 = beta0/180.0*3.14;
    //====== obtain the orinatation of the beta
    Vec3D P1 (v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos()) ;
    Vec3D P2 (v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos()) ;
    double l = (P2-P1).norm();
    Vec3D l1 = N2-N1+(P2-P1)*(1.0/l);
    if(l1.norm()<1)
        beta=-beta;
    E = -e0-e1*cos(N*(m_Angle2D-theta0))+e2*(beta-beta0)*(beta-beta0);

    return E;
}
double Energy::F10(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = -e0+e1*cosN(phi-phi0)+e2*exp(-alpha(theta-theta0))
    double E = 0;
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double phi0 = var.at(3);
    double e2 = var.at(4);
    double alpha = var.at(5);
    double theta0 = var.at(6);
    
    
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";

    Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D geodesic_dir=(X2-X1);
    Vec3D *pBox=v1->GetBox();
    
    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
        }
    }
    
    
    Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
    y1(2)=0;
    y1=y1*(1/(y1.norm()));
    Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
    y2(2)=0;
    y2=y2*(1/(y2.norm()));
    Vec3D n(0,0,1);
    Vec3D d1 = (v1->GetInclusion())->GetLDirection();
    Vec3D d2 = (v2->GetInclusion())->GetLDirection();
    double cos1 = y1.dot(y1,d1);
    double sin1 = n.dot(n*y1,d1);
    double cos2 = y1.dot(y2,d2);
    double sin2 = n.dot(n*y2,d2);
    double S_an = sin1*sin2+cos1*cos2;
    m_Angle2D=acos(S_an);
    
    
    Vec3D gd1 = (v2->GetL2GTransferMatrix())*d1;
    Vec3D gd2 = (v2->GetL2GTransferMatrix())*d2;
    
    m_Angle2D  = acos(n.dot(gd1,gd2));
    theta0 = theta0/180.0*3.14;
    
    E = e0+e1*cos(N*(m_Angle2D-phi0))+e2*exp(-alpha*(m_Angle2D-theta0)*(m_Angle2D-theta0));
    E = -E;
    
    return E;
}
double Energy::F11(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = e0+e1*cosN(phi-phi0)+e2*(3*(m1.r)(m2.r2)-m1.m2)
    
    /// m1 = Dir1+DeltaD
    double E = 0;
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double phi0 = var.at(3);
    double e2 = var.at(4);
    double Q0 = var.at(5);
    

    
    
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";
    
    Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D geodesic_dir=(X2-X1);
    Vec3D *pBox=v1->GetBox();
    
    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
        }
    }
    double GEOLength = geodesic_dir.norm();
    Vec3D GEOUnit = geodesic_dir*(1.0/GEOLength);   /// geo direction
    
    Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
    y1(2)=0;
    y1=y1*(1/(y1.norm()));
    Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
    y2(2)=0;
    y2=y2*(1/(y2.norm()));
    Vec3D n(0,0,1);
    Vec3D d1 = (v1->GetInclusion())->GetLDirection();
    Vec3D d2 = (v2->GetInclusion())->GetLDirection();
    double cos1 = y1.dot(y1,d1);
    double sin1 = n.dot(n*y1,d1);
    double cos2 = y1.dot(y2,d2);
    double sin2 = n.dot(n*y2,d2);
    double S_an = sin1*sin2+cos1*cos2;
    m_Angle2D=acos(S_an);
    
    Vec3D q0(0,0,tan(Q0));
    Vec3D gd1 = d1+q0;
    Vec3D gd2 = d2+q0;
    gd1 = gd1*(1/(gd1.norm()));
    gd2 = gd2*(1/(gd2.norm()));
    gd1 = (v1->GetL2GTransferMatrix())*gd1;
    gd2 = (v2->GetL2GTransferMatrix())*gd2;
    
    double Q1=gd1.dot(gd1,GEOUnit);
    double Q2=gd2.dot(gd2,GEOUnit);
    
    
    E = e0+e1*cos(N*(m_Angle2D-phi0))+e2*(3*Q1*Q2-n.dot(gd1,gd2));
    E = -E;
    
    return E;
}
double Energy::TwoVectorFieldInteractionEnergy(int vf_layer, links * p_edge) {

    vertex * p_v1 = p_edge->GetV1();
    vertex * p_v2 = p_edge->GetV2();
    
    if(vf_layer<0 || vf_layer >= p_v1->GetNumberOfVF()){
        std::cout<<"---> error: unexpected, this function should not be called \n";
        return 0;
    }
//--> now use the type if to get the interaction parameters of this pars
    double e_int = 0;
    VectorField *VF1 =  p_v1->GetVectorField(vf_layer);
    VectorField *VF2 =  p_v2->GetVectorField(vf_layer);
    Vec3D d1 = VF1->GetLDirection();
    Vec3D d2 = VF2->GetLDirection();
    int id1 = VF1->GetInclusionType()->ITid;
    int id2 = VF2->GetInclusionType()->ITid;
    PairInt pair_ab = m_pInt->GetPairInt(id1,id2);
    std::vector <double> ff = pair_ab.Varibale;
    //---> get the type of function that they interact, it is an intger
    int FunctionType  = pair_ab.FunctionType;
   // std::cout<<ff.size()<<"  ff.size \n";
   // m_pInt->Print();

    switch (FunctionType) {
        case 0:{
            e_int = 0.0;
             break;
        }
         case 1: {
             
             double theta = (ff[2] != 0) ? AngleDiff_ParallelTransport(d1, d2, p_edge) : 0.0;
             e_int = -ff[1] + ff[2] * cos(double(ff[0]) * theta);
             break;
         }
        default:{
             std::cerr << "---> Error: Unrecognized function type ID for vector field--> " << FunctionType << std::endl;
             exit(0);
        }
     }
    
    if(p_edge->GetMirrorFlag()) {
        p_edge->UpdateVFIntEnergy(vf_layer, e_int/2.0);
        (p_edge->GetMirrorLink())->UpdateVFIntEnergy(vf_layer, e_int/2.0);
    }
    else {
        // divided by 2 here is also fine as everywhere will be muliplied by two again
        p_edge->UpdateVFIntEnergy(vf_layer, e_int/2.0);
    }
    return e_int;
}
double Energy::AngleDiff_ParallelTransport(Vec3D &d1, Vec3D &d2, links* p_link) {
    /*
     * @brief Calculate the angle between two vectors after parallel transport along the geodesic direction.
     *
     * This function computes the angle between two vectors after parallel transport along the geodesic direction
     * connecting two vertices. The geodesic direction is determined by subtracting the position vectors of the
     * two vertices. The resulting angle represents the deviation between the directions of the two vectors
     * after transporting them along the geodesic path.
     */
    
        vertex *v1 = p_link->GetV1();
        vertex *v2 = p_link->GetV2();

        Vec3D X1 = v1->GetPos();
        Vec3D X2 = v2->GetPos();
        Vec3D geodesic_dir=(X2-X1);
    
        for (int i=0;i<3;i++)
        {
            if(fabs(geodesic_dir(i))>  m_Box(i)/2)
            {
                if(geodesic_dir(i)<0)
                    geodesic_dir(i) = geodesic_dir(i) + m_Box(i);
                else if(geodesic_dir(i)>0)
                    geodesic_dir(i) = geodesic_dir(i) - m_Box(i);
            }
        }
        Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
        y1(2)=0;
        y1.normalize();
        Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
        y2(2)=0;
        y2.normalize();
        Vec3D n(0,0,1);
        double cos1 = y1.dot(y1,d1);
        double sin1 = n.dot(n*y1,d1);
        double cos2 = y1.dot(y2,d2);
        double sin2 = n.dot(n*y2,d2);
    
        // Compute theta
        double S_an = sin1 * sin2 + cos1 * cos2;
        double theta = acos(S_an);

    return theta;
}
double Energy::Geo_Theta(vertex *v1, vertex *v2) {
    /*
     * @brief Calculate the angle between two vectors after parallel transport along the geodesic direction.
     *
     * This function computes the angle between two vectors after parallel transport along the geodesic direction
     * connecting two vertices. The geodesic direction is determined by subtracting the position vectors of the
     * two vertices. The resulting angle represents the deviation between the directions of the two vectors
     * after transporting them along the geodesic path.
     *
     * @param v1 Pointer to the first vertex.
     * @param v2 Pointer to the second vertex.
     * @return The angle (in radians) between the transported vectors.
     */
    
        Vec3D X1 = v1->GetPos();
        Vec3D X2 = v2->GetPos();
        Vec3D geodesic_dir=(X2-X1);
    
        for (int i=0;i<3;i++)
        {
            if(fabs(geodesic_dir(i))>  m_Box(i)/2)
            {
                if(geodesic_dir(i)<0)
                    geodesic_dir(i) = geodesic_dir(i) + m_Box(i);
                else if(geodesic_dir(i)>0)
                    geodesic_dir(i) = geodesic_dir(i) - m_Box(i);
            }
        }
        Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
        y1(2)=0;
        y1.normalize();
        Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
        y2(2)=0;
        y2.normalize();
        Vec3D n(0,0,1);
        Vec3D d1 = (v1->GetInclusion())->GetLDirection();
        Vec3D d2 = (v2->GetInclusion())->GetLDirection();
        double cos1 = y1.dot(y1,d1);
        double sin1 = n.dot(n*y1,d1);
        double cos2 = y1.dot(y2,d2);
        double sin2 = n.dot(n*y2,d2);
    
        // Compute theta
        double S_an = sin1 * sin2 + cos1 * cos2;
        double theta = acos(S_an);

    return theta;
}
double Energy::CalculateVectorFieldMembraneBindingEnergy(VectorField* p_vf, vertex *p_vertex){
    
    double en = 0;
    int layer = p_vf->GetLayer();
    en += m_pState->GetExternalFieldOnVectorFields()->GetCouplingEnergy(layer, p_vf, p_vertex);

    Vec3D l_direction = p_vf->GetLDirection();
    InclusionType* inc_type = p_vf->GetInclusionType();
    
if(p_vertex->GetVertexType() == 0) {
    double c1 = p_vertex->GetP1Curvature();
    double c2 = p_vertex->GetP2Curvature();

    double mean_times2 = c1 + c2;
    double gussian = c1 * c2;
    double area = p_vertex->GetArea();
    
    double k0 = inc_type->ITk;
    double kg = inc_type->ITkg;
    double k1 = inc_type->ITk1;
    double k2 = inc_type->ITk2;
    double c0 = inc_type->ITc0;
    double cp10 = inc_type->ITc1;
    double cn20 = inc_type->ITc2;
    //--- kappa/2*(2H-c0)^2-kgK+
    double ev = k0 * (mean_times2-c0) * (mean_times2-c0) - kg * gussian;
    en += ev * area;
    
//--> if the k1 and k2 are zero, we do not need to calculate the rest
    if(k1 != 0 || k2 != 0){  // k1/2(cp-cp0)^2+k2/2(cn-cn0)^2
    
        double Cos = l_direction(0);
        double Sin = l_direction(1);
        double Cp = c1 * Cos * Cos + c2 * Sin * Sin;
        double Cn = c2 * Cos * Cos + c1 * Sin * Sin;
        double Delta_Cp = Cp - cp10;
        double Delta_Cn = Cn - cn20;
        en += (k1 * Delta_Cp * Delta_Cp + k2 * Delta_Cn * Delta_Cn) * area;
    }
}
else if( p_vertex->GetVertexType() == 1) {
    
    double geo_c = p_vertex->GetGeodesicCurvature();
    double norm_c = p_vertex->GetNormalCurvature();
    double length = p_vertex->GetLength();
    
        double lambda = inc_type->ITelambda;
        double kg     = inc_type->ITekg;
        double kn     = inc_type->ITekn;
        double cn0    = inc_type->ITecn;
        
        en += (lambda + kg*geo_c*geo_c) * length;
        
        if(kn != 0) {
            
            double Cos = l_direction(0);
            kn = kn*Cos*Cos;
            double Delta_norm_c = (norm_c-cn0);
            en += (kn*Delta_norm_c*Delta_norm_c)*length;
        }
    
}
else{
    
    std::cout<<" error-> this is unexpected 903 \n";
}
    p_vf->UpdateMembraneBindingEnergy(en);
    
    return en;
}
double Energy::CalculateVectorFieldMembraneBindingEnergy(vertex *p_vertex){
    
    double T_en = 0;
    int no_vf = p_vertex->GetNumberOfVF();
    
    for(int i = 0; i<no_vf; i++) {
        T_en += CalculateVectorFieldMembraneBindingEnergy(p_vertex->GetVectorField(i), p_vertex);
    }
    
    return T_en;
    
}
std::string Energy::CurrentState(){
    
    std::string state = AbstractEnergy::GetBaseDefaultReadName() + " = " + GetDerivedDefaultReadName();
    state = state + "\n Kappa = "+Nfunction::D2S(2*m_kappa)+" "+Nfunction::D2S(m_kappa_G)+" "+Nfunction::D2S(m_SCurvature0);
    state = state + "\n Edge_Parameters = "+Nfunction::D2S(m_Lambda)+" "+Nfunction::D2S(m_Kappa_Geo)+" "+Nfunction::D2S(m_Kappa_Norm);
    state = state + "\n VertexArea = "+Nfunction::D2S(m_Ka)+" "+Nfunction::D2S(m_Area0/sqrt(3)-0.5)+" "+Nfunction::D2S(m_Kl)+" "+Nfunction::D2S((m_l0*m_l0-1)/2);

    return state;
}
