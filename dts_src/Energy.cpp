

#include <time.h>
#include "Energy.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
Energy::Energy(Inclusion_Interaction_Map * pint)
{
m_pInt=pint;
}
Energy::~Energy()
{
    
}
double Energy::SingleVertexEnergy(vertex *pv)
{

    std::vector<double> Curve=pv->GetCurvature();
    double kappa=pv->GetKappa();     //   this kappa is kappa/2  we have done it simulation class
    double kappag=pv->GetKappaG();     //
    double area=pv->GetArea();
    double mean=(Curve.at(0)+Curve.at(1));
    double gussian=(Curve.at(0))*(Curve.at(1));
    double Energy=0.0;

    if(pv->VertexOwnInclusion()==true)
    {
        inclusion *inc=pv->GetInclusion();
        InclusionType *inctype = inc->GetInclusionType();
        double k0 = inctype->ITk;
        double kg = inctype->ITkg;
        double k1 = inctype->ITk1;
        double k2 = inctype->ITk2;
        double c0 = inctype->ITc0;
        double c1 = inctype->ITc1;
        double c2 = inctype->ITc2;

            double H=(mean-c0);
            k0=k0/2.0;
            Energy+=(k0*H*H-kg*gussian)*area;         /// this means that inclsuion overwrite the vertex bending rigidity
        if(k1!=0 || k2!=0)
        {
            Vec3D LD=inc->GetLDirection();
            k1=k1/2.0;
            k2=k2/2.0;
            double Cos=LD(0);
            double Sin=LD(1);
            double C1=Curve.at(0)*Cos*Cos+Curve.at(1)*Sin*Sin;
            double C2=Curve.at(1)*Cos*Cos+Curve.at(0)*Sin*Sin;
            double H1=(C1-c1);
            double H2=(C2-c2);
            Energy+=(k1*H1*H1+k2*H2*H2)*area;

        }
    }
    else
    {
        Energy=(kappa*mean*mean-kappag*gussian)*area;

    }
    

    pv->UpdateEnergy(Energy);


#if TEST_MODE == Enabled
    if(isnan(Energy)==1)
    std::cout<<Energy<<" energy not got \n";
#endif
    return Energy;
}
double Energy::TotalEnergy(std::vector<vertex *> pVeretx, std::vector<links *> pLinks)
{
    double E=0.0;
    
    for (std::vector<vertex *>::iterator it = pVeretx.begin() ; it != pVeretx.end(); ++it)
    {
        E+=SingleVertexEnergy(*it);
    }
    for (std::vector<links *>::iterator it = pLinks.begin() ; it != pLinks.end(); ++it)
    {
        E+=TwoInclusionsInteractionEnergy(*it);
    }
    
    return E;
}
double Energy::Energy_OneLinkFlip(links * plinks)
{
 double E=0.0;
    
    links *mirorl=plinks->GetMirrorLink();
    vertex *v1=plinks->GetV1();
    vertex *v2=plinks->GetV2();
    vertex *v3=plinks->GetV3();
    vertex *v4=mirorl->GetV3();

    
    
    E+=SingleVertexEnergy(v1);
    E+=SingleVertexEnergy(v2);
    E+=SingleVertexEnergy(v3);
    E+=SingleVertexEnergy(v4);
    
    
    

    
    
    
  return E;
    

}
double Energy::Energy_OneVertexMove(vertex * pVeretx)
{

    double E=0.0;
    std::vector<vertex *> NpVer=pVeretx->GetVNeighbourVertex(); /// Get The vertexs on the ring


    
    
    E+=SingleVertexEnergy(pVeretx);

    for (std::vector<vertex *>::iterator it = NpVer.begin() ; it != NpVer.end(); ++it)
    {
        E+=SingleVertexEnergy(*it);
    }
    

    return E;
}
double Energy::TwoInclusionsInteractionEnergy(links * lp)
{

    vertex * v1 = lp->GetV1();
    vertex * v2 = lp->GetV2();
    double e=0;
    
    bool has1=v1->VertexOwnInclusion();
    bool has2=v2->VertexOwnInclusion();
    
    if(has1==true && has2==true)
    {
        int id1=(v1->GetInclusion())->GetInclusionTypeID();
        int id2=(v2->GetInclusion())->GetInclusionTypeID();
        PairInt pair_ab = m_pInt->GetPairInt(id1,id2);
        std::vector <double> ff = pair_ab.Varibale;
        int FunctionType  = pair_ab.FunctionType;
        
            double theta = 0;
         if(FunctionType == 0)
         {
                e=0;
         }
          else if(FunctionType == 1)
          {
            if(ff.at(2)!=0)
            theta = Geo_Theta(v1,v2);
            e=InteractionFunction(ff.at(0), ff.at(1),ff.at(2),theta);
          }
         else if(FunctionType == 10)
         {
             m_Angle3D = 0;
             m_Angle2D = 0;
             e = F10(v1,v2,ff);
         }
        else if(FunctionType == 11)
        {
            m_Angle3D = 0;
            m_Angle2D = 0;
            e = F11(v1,v2,ff);
        }
         else
         {
             std::cout<<"---> Error: Unregognized function typeid -->"<<FunctionType<<std::endl;
             exit(0);
         }

    }
    
    
    lp->UpdateIntEnergy(e/2.0);
    if(lp->GetMirrorFlag()==true)      //// this check is not needed for now, for later developments 
    (lp->GetMirrorLink())->UpdateIntEnergy(e/2.0);
    

    
    return e;
}
double Energy::InteractionFunction(double N, double A, double B, double theta)
{
    double e = 0;


        e=1+cos(double(N)*theta);
        e=-A+B*e;

    return e;
}
double Energy::Geo_Theta(vertex *v1, vertex *v2)
{

    double theta;
        
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
         theta=acos(S_an);
        
    
    
    return theta;
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
    
    
    /*
      /
     /Q0
     -------->
     */
    
    
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
