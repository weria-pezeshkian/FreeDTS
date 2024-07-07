#include "Constant_NematicForceByVectorFields.h"


Constant_NematicForceByVectorFields::Constant_NematicForceByVectorFields(std::string data_stream) {
    
    std::vector<std::string> data = Nfunction::Split(data_stream);
    if(data.size() == 0){
        std::cout<<"---> error: data provided for constant nematic force is not enough \n";
    }
    for (std::vector<std::string>::iterator it = data.begin() ; it != data.end(); ++it){
        double f = Nfunction::String_to_Double(*it);
        m_F0.push_back(f);
    }
    m_ActiveEnergy = 0;
}
Constant_NematicForceByVectorFields::~Constant_NematicForceByVectorFields() {
    
}
double Constant_NematicForceByVectorFields::Energy_of_Force(vertex *p_vertex, Vec3D dx) {
    
    if( p_vertex->GetNumberOfVF() == 0){
        return 0;
    }
    
    if(  m_F0.size() != p_vertex->GetNumberOfVF() ){
        std::cout<<"---> error: data provided for constant nematic force is not enough \n";
        return 0;
    }
    double En = 0;
    for (int i = 0; i< p_vertex->GetNumberOfVF(); i++){
        
        if(m_F0[i] !=0 ){
            std::vector <links *> nl = p_vertex->GetVLinkList();
            Vec3D Force;
            for (std::vector<links *>::iterator it = nl.begin() ; it != nl.end(); ++it) {
                vertex *p_vertex2 = (*it)->GetV2();
                Force = Force + ActiveNematicForce_1(i, p_vertex2, p_vertex);
            }
       
            Tensor2  G2L = p_vertex->GetG2LTransferMatrix();
            Vec3D ldx = G2L * dx;
            En += m_F0[i] * Vec3D::dot(ldx , Force);  // in the local space F = -zeta*div(Q); E=zeta*div(Q)*dX
        }
    }
    return En;
}
Vec3D Constant_NematicForceByVectorFields::ActiveNematicForce_1(int layer, vertex *v2, vertex *v1) // gives force in the local coordinate
{
    Vec3D f;
    
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

        Vec3D nem1_d = (v1->GetVectorField(layer))->GetLDirection();
        Vec3D nem2_d = (v2->GetVectorField(layer))->GetLDirection(); // this needs to be transported to v1
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
std::string Constant_NematicForceByVectorFields::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    for (std::vector<double>::iterator it = m_F0.begin() ; it != m_F0.end(); ++it){
        state += "  "+ Nfunction::D2S(*it);
    }
    return state;
}

