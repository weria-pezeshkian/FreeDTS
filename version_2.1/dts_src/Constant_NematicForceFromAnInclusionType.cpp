#include "Constant_NematicForceFromAnInclusionType.h"


Constant_NematicForceFromAnInclusionType::Constant_NematicForceFromAnInclusionType(double f0, std::string inc_typename) : m_IncType(inc_typename){
    m_F0 = f0;
    m_ActiveEnergy = 0;



}
Constant_NematicForceFromAnInclusionType::~Constant_NematicForceFromAnInclusionType() {
    
}
double Constant_NematicForceFromAnInclusionType::Energy_of_Force(vertex *pv, Vec3D dx) {
    
    if(!pv->VertexOwnInclusion() ||  m_F0 == 0 )
        return 0;
    
    if(pv->GetInclusion()->GetInclusionType()->ITName != m_IncType )
        return 0;
    
    double En = 0;
    
        std::vector <links *> nl = pv->GetVLinkList();
        Vec3D Force;
        for (std::vector<links *>::iterator it = nl.begin() ; it != nl.end(); ++it)
        {
            vertex *pv2 = (*it)->GetV2();
            if(pv2->VertexOwnInclusion()){
                Force = Force + ActiveNematicForce_1(pv2, pv);
            }
        }
        Tensor2  G2L = pv->GetG2LTransferMatrix();
        Vec3D ldx = G2L * dx;
        En = m_F0 * Vec3D::dot(ldx , Force);  // in the local space F = -zeta*div(Q); E=zeta*div(Q)*dX
    
    return En;
}
Vec3D Constant_NematicForceFromAnInclusionType::Inclusion_Force(vertex *pv) {
    
    Vec3D TotalForce(0,0,0);
    if(!pv->VertexOwnInclusion() ||  m_F0 == 0 )
        return TotalForce;
    
    
        std::vector <links *> nl = pv->GetVLinkList();
        for (std::vector<links *>::iterator it = nl.begin() ; it != nl.end(); ++it)
        {
            vertex *pv2 = (*it)->GetV2();
            if(pv2->VertexOwnInclusion()){
                TotalForce = TotalForce + ActiveNematicForce_1(pv2, pv);
            }
        }
        Tensor2  L2G = pv->GetL2GTransferMatrix();
        TotalForce = (L2G * TotalForce) * m_F0;  // in the local space F = -zeta*div(Q); E=zeta*div(Q)*dX
    
    return TotalForce;
}
Vec3D Constant_NematicForceFromAnInclusionType::ActiveNematicForce_1(vertex *v2, vertex *v1) // gives force in the local coordinate
{
    Vec3D f;
    
    if(!v2->VertexOwnInclusion())
        return f;
    
        // obtaining geodesic direction
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
        double d = geodesic_dir.norm(); // obtaining the distance between the two vertex
        Vec3D gdir=(v1->GetG2LTransferMatrix())*geodesic_dir;
        gdir(2)=0;
        gdir = gdir * (1/(gdir.norm())); // unit vector d in the theory
        // end obtaining geodesic direction
        Vec3D gdir_2 = (v2->GetG2LTransferMatrix())*geodesic_dir;
        gdir_2(2)=0;
        gdir_2 = gdir_2 * (1/(gdir_2.norm())); // unit vector d in the theory
        Vec3D n(0,0,1);

        Vec3D nem1_d = (v1->GetInclusion())->GetLDirection();
        Vec3D nem2_d = (v2->GetInclusion())->GetLDirection(); // this needs to be transported to v1
        // claculate div Q contribution
     
        // we only need sin and cosine
        double cos1 = nem1_d.dot(nem1_d,gdir);
        double sin1 = n.dot(n*gdir,nem1_d);
        double cos2 = nem2_d.dot(nem2_d,gdir_2);
        double sin2 = n.dot(n * gdir_2 , nem2_d);
        double sinDeltaT = sin2*sin1+cos2*cos1;  //cos (theta2-theta1) after parallel transport
        double Delta_T = acos(sinDeltaT);
        Delta_T = -2*Delta_T/d;
        
        double sin2T = sin1*cos1+cos1*sin1;  //sin (2theta1)
        double cos2T = cos1*cos1-sin1*sin1;  //cos (2theta1)
        double cosphi = gdir(0);
        double sinphi = gdir(1);



        f(0) = sin2T * cosphi - cos2T * sinphi;   // sin(2theta-phi)
        f(1) = -sin2T * sinphi - cos2T * cosphi;   // -cos(2theta-phi)

        f = f*Delta_T;
    

    return f;
}
std::string Constant_NematicForceFromAnInclusionType::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_F0);
    return state;
}

