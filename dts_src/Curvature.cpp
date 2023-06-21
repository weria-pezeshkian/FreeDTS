


#include "Curvature.h"
#include "Tensor2.h"
Curvature::Curvature()
{
    
}
Curvature::~Curvature()
{

}
void Curvature::CalculateCurvature(vertex * pvertex)
{
    m_pVertex=pvertex;
    std::vector<triangle *> Ntr=m_pVertex->GetVTraingleList();
    
    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        Vec3D v=(*it)->GetNormalVector();
        double a=(*it)->GetArea();
        v=v*a;
        Normal=Normal+v;
        Area+=a;
    }

    Area=Area/3.0;
   

  
    if(Area==0)
    {
    std::cout<<Ntr.size()<<"\n";
	std::string sms=" error----> vertex has a zero area \n";
	std::cout<<sms<<"\n";
        exit(0);
    }
    else if(Area<0)
    {
	std::string sms=" error----> vertex has a negetive area \n";
	std::cout<<sms<<"\n";
        exit(0);
    }


    double no=Normal.norm();
    no=1.0/no;
    Normal=Normal*no;
    m_pVertex->UpdateNormal_Area(Normal,Area);
        ///=======
    //=== Shape Operator
    //========
    Tensor2  SV;
    Tensor2 IT('I');
    
    Tensor2 P=IT-IT.makeTen(Normal);

    Tensor2 Pt=P.Transpose(P);
    std::vector<links *> NLinks=m_pVertex->GetVLinkList();

    for (std::vector<links *>::iterator it = NLinks.begin() ; it != NLinks.end(); ++it)
    {
       if((*it)->GetMirrorFlag()==true)
	   {
        Vec3D ve=(*it)->GetNormal();
        double we=ve.dot(Normal,ve);
        

        Vec3D Be=(*it)->GetBe();
        double he=(*it)->GetHe();
        Vec3D Se;
        for (int n=0;n<3;n++)
        {
            Se(0)=Se(0)+P(0,n)*Be(n);
            Se(1)=Se(1)+P(1,n)*Be(n);
            Se(2)=Se(2)+P(2,n)*Be(n);
        }
        
        double ff=Se.norm();
        
        
        if(ff==0)
        {
            std::cout<<"-----> Error: projection is zero error"<<"\n";
            exit(0);
        }
        else
        {
            Se=Se*(1.0/ff);
        }
        
        Tensor2 Q=P.makeTen(Se);
        
        SV=SV+(Q)*(we*he);
	}
    }
    ///=============
    //==== Find Curvature and local frame
    //=============
    
    Tensor2 Hous = Householder(Normal);
    Tensor2 LSV;// Local SV
    LSV=(Hous.Transpose(Hous))*(SV*Hous);    // LSV is the curvature matrix in the local frame, the it is a 2*2 minor matrix since all the 3 component are zero.

    double b=LSV(0,0)+LSV(1,1);
    double c=LSV(0,0)*LSV(1,1)-LSV(1,0)*LSV(0,1);


    double delta=b*b-4*c;
    double c1,c2;
    if(delta>0.0)
    {
    delta=sqrt(delta);
    c1=b+delta;
    c1=0.5*c1;      // c1 always will be larger then c2
    c2=b-delta;
    c2=0.5*c2;
        

    }
    else if (fabs(delta)<0.0001)
    {
        
    c1=0.5*b;
    c2=c1;

    }
    else
    {
        c1=100;
        c2=100;
        std::cout<<"WARNING: faild to find curvature on vertex "<<pvertex->GetVID()<<"  because delta is "<<delta<<"  c1 and c2 are set to 100 \n";
        std::cout<<" if you face this too much, you should stop the job and .... \n";
    }

    

   // if(true) in general we do not need this, only if we have directional inclsuions
    {
    Tensor2 EigenvMat('O');
    
    double p=LSV(0,0);
    double q=LSV(0,1);


    double size=sqrt(q*q+(c1-p)*(c1-p));                   // The Eigenvectors can be calculated using this equation LSV*R=c1*R
    EigenvMat(0,0)=q/size;                                  // only one of them needs to be calculated, one is normal vector and the other is perpendicular to first one
    EigenvMat(1,0)=(c1-p)/size;
    EigenvMat(0,1)=-EigenvMat(1,0);
    EigenvMat(1,1)=EigenvMat(0,0);
    EigenvMat(2,2)=1;
        


        // std::cout<<(Hous*Normal)(0)<<"   "<<(Hous*Normal)(1)<<"   "<<(Hous*Normal)(2)<<"   \n";

 // Tensor2 TransferMatGL=EigenvMat*Hous;   /// This matrix transfers vectors from Global coordinate to local coordinate//
 // Tensor2 TransferMatLG = TransferMatGL.Transpose(TransferMatGL);//// This matrix transfers vectors from local coordinate to global coordinate
        
        
        ///  this is correct, We can check by applying transpose(E)*t1 = (1,0,0)
        
     Tensor2 TransferMatLG=Hous*EigenvMat;   /// This matrix transfers vectors from local coordinate to global coordinate
        

        
     Tensor2 TransferMatGL=TransferMatLG.Transpose(TransferMatLG);   /// This matrix transfers vectors from Global coordinate to local coordinate
        
        m_pVertex->UpdateL2GTransferMatrix(TransferMatLG);
        m_pVertex->UpdateG2LTransferMatrix(TransferMatGL);
        
        
        
        if(pvertex->VertexOwnInclusion()==true)
        {
            inclusion * in=pvertex->GetInclusion();
            Vec3D LD=in->GetLDirection();
        }
        

//
 #if TEST_MODE == Enabled
      Tensor2 Before = TransferMatLG;
        Vec3D P1(1,0,0);
        Vec3D P2(0,1,0);
        
        P1=TransferMatLG*P1;
        P2=TransferMatLG*P2;

        Vec3D p1=P1;
        if(P1.dot(P1*P2,Normal)<0)
        P1=P1*(-1);
        
        Tensor2 Q(P1,P2,Normal);
        
        TransferMatGL = Q;
        
        
        TransferMatLG = Q.Transpose(TransferMatGL);
        
        Tensor2 After = TransferMatLG;

  
        
        Vec3D a1(1,0,0);
        Vec3D a2(0,1,0);
        Vec3D a3(0,0,1);

        if(P1.dot(p1*P2,Normal)<0)
        {
            std::cout<<((TransferMatLG*a1)*(TransferMatLG*a2))(0)<<"   "<<((TransferMatLG*a1)*(TransferMatLG*a2))(1)<<"   "<<((TransferMatLG*a1)*(TransferMatLG*a2))(2)<<"\n";
            
           std::cout<<(TransferMatGL*(TransferMatLG*a3))(0)<<"   "<<(TransferMatGL*(TransferMatLG*a3))(1)<<"   "<<(TransferMatGL*(TransferMatLG*a3))(2)<<"\n";

            std::cout<<"=================\n";
            std::cout<<Before(0,0)<<"   "<<Before(0,1)<<"   "<<Before(0,2)<<"\n";
            std::cout<<Before(1,0)<<"   "<<Before(1,1)<<"   "<<Before(1,2)<<"\n";
            std::cout<<Before(2,0)<<"   "<<Before(2,1)<<"   "<<Before(2,2)<<"\n";
            
            std::cout<<"====== after ===========\n";
            std::cout<<After(0,0)<<"   "<<After(0,1)<<"   "<<After(0,2)<<"\n";
            std::cout<<After(1,0)<<"   "<<After(1,1)<<"   "<<After(1,2)<<"\n";
            std::cout<<After(2,0)<<"   "<<After(2,1)<<"   "<<After(2,2)<<"\n";
            
            std::cout<<"\n";
            

        }
        ///// here should be corrected  . we need a good test to see if the transformation is correct
        Vec3D localN(0,0,1);
        Vec3D GN1 = TransferMatLG*localN;
    if( fabs(Normal(2)-GN1(2))>0.001 || fabs(Normal(1)-GN1(1))>0.001 || fabs(Normal(0)-GN1(0))>0.001)
    {
        std::cout<<" error code 1110319: something is unusuall \n";
        std::cout<<" transfer matrix does not work well \n";
        std::cout<<Normal(0)<<"   "<<Normal(1)<<"   "<<Normal(2)<<"   G \n";
        std::cout<<GN1(0)<<"   "<<GN1(1)<<"   "<<GN1(2)<<"  Local to G  \n";

        std::cout<<"\n ";
    }
        Vec3D LN = TransferMatGL*Normal;
        if( fabs(LN(0))>0.0001 || fabs(LN(1))>0.0001 || fabs(1-LN(2))>0.001)
        {
            std::cout<<" error code 1110319: something is unusuall \n";
            std::cout<<" transfer matrix does not work well \n";
            std::cout<<LN(0)<<"   "<<LN(1)<<"   "<<LN(2)<<"  Local to G  \n";
            
            std::cout<<"\n ";
        }

        
         Vec3D k1(1,2,3);
         Vec3D k2=TransferMatLG*k1;
         Vec3D k3=TransferMatGL*k2;
        
        if(fabs(k1(0)-k3(0))>0.0001 ||fabs(k1(1)-k3(1))>0.0001  || fabs(k1(2)-k3(2))>0.0001 )
        {
            std::cout<<" transfer matrix does not work well: when we transfer and back v=1,2,3 \n";
         std::cout<<"===============\n";
         std::cout<<k1(0)<<" "<<k1(1)<<" "<<k1(2)<<" \n";
         std::cout<<k3(0)<<" "<<k3(1)<<" "<<k3(2)<<" \n";
        }
#endif
        
        
    }


    
    c1=c1/Area;
    c2=c2/Area;
    m_pVertex->UpdateCurvature(c1,c2);
    
    

 
}
Tensor2 Curvature::Householder(Vec3D N)
{
    
    Tensor2 Hous;
    Vec3D Zk;
    Zk(2)=1.0;
    Zk=Zk+N;
    Zk=Zk*(1.0/Zk.norm());
    
    Tensor2 I('I');
    Tensor2 W=Hous.makeTen(Zk);
    Hous=(I-W*2)*(-1);
    
    return Hous;

    

   
    /*
    Zk(0) = 0;Zk(1) = 0;
    Zk(2)=1.0;
    double SignT=1;
    if((1+N(2))>(1-N(2)))   // 1+/-Normal(2) is always larger then 1
    {
     Zk=Zk+N;
    SignT=-1;
     }
    else if((1+N(2))<=(1-N(2)))
     Zk=Zk-N;
    Zk=Zk*(1.0/Zk.norm());
    W=Hous.makeTen(Zk);
    Hous=(I-W*2)*(SignT);
    
    
    std::cout<<"====== old ===========\n";
    std::cout<<Hous(0,0)<<"   "<<Hous(0,1)<<"   "<<Hous(0,2)<<"\n";
    std::cout<<Hous(1,0)<<"   "<<Hous(1,1)<<"   "<<Hous(1,2)<<"\n";
    std::cout<<Hous(2,0)<<"   "<<Hous(2,1)<<"   "<<Hous(2,2)<<"\n";*/
   // double SignT=1;
   
}
/// normal vector update






