

#include <time.h>
#include "src/Analysis.h"
#include "src/LMatrix.h"
#include "src/WritevtuFiles.h"
#include "src/Nfunction.h"

struct L_CLUSTER {
    
    std::vector<vertex *> Vertex;
    std::vector<int> NCluster;
    std::vector<vertex *> ItsV;
    std::vector<vertex *> NV;
    int checked;
    Vec3D LX;
} ;

Analysis::Analysis(State *pState)
{
    Nfunction f;
    m_TRJ = pState->GetTrajectory();
    std::cout<<m_TRJ.size()<<"size of trj \n";
    m_pState =pState;
    double k= pState->m_BendingRigidity;
    double kg= pState->m_GaussianRigidity;
    int i=0;
    m_totEn = 0;
    m_totVolume = 0;
    m_totArea = 0;
    std::ofstream reponsedata;
    reponsedata.open("response_data.xvg");
    reponsedata<<"@@    step   E      V        totArea   A_out A_in     \n";

    for (std::vector<MeshBluePrint>::iterator it = m_TRJ.begin() ; it != m_TRJ.end(); ++it)
    {
        MESH mesh;
        mesh.GenerateMesh(*it, k, kg);

        UpdateEnergy(mesh);
        UpdateVolume(mesh);
       // RemovePBC(mesh);
        //Center(mesh);
        SeparateInOutByVer(mesh);
        std::vector<std::vector<double> > data;
        std::vector<std::string> num_name;
        // include data
        std::vector<vertex*> pV = mesh.m_pAllV;
     //   std::vector<trinagle*> pT = mesh.m_Triangle;

        //===============================
        
        
        
        std::vector<double> Surf; // mean curvature
        std::vector<double> en; // energy

        std::string name1="Surf";
        num_name.push_back(name1);
        std::string name2="En";
        num_name.push_back(name2);

        for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
        {
            std::vector<double> c = (*it)->GetCurvature();
            en.push_back((*it)->GetEnergy());
            Surf.push_back((*it)->GetGroup());

        }
        data.push_back(Surf);
        data.push_back(en);



        double a_in = 0;
        double a_out = 0;
      
        for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
        {
		if((*it)->GetGroup()==0)
		a_out+= (*it)->GetArea();
		a_in+= (*it)->GetArea();

        }

        //=====================================
        
        reponsedata<<i<<"   "<<m_totEn<<"      "<<m_totVolume<<"      "<<m_totArea<<"  "<<a_out<<"  "<<a_in<<"   \n";
        
        
        //======== writing vtu files
        WritevtuFiles vtu(m_pState);
        std::string fi="conf"+f.Int_to_String(i)+".vtu";
        vtu.Writevtu_Plus_NumberList(num_name, data, mesh.m_pAllV,mesh.m_pAllT,mesh.m_pHL,fi);
        i++;
    }

    

    
    
    
    
}
Analysis::~Analysis()
{
}
void Analysis::UpdateVolume(MESH &mesh)
{
    std::vector<triangle *> pT=mesh.m_pAllT;
    double V=0.0;
    double A=0.0;

    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
    {
        vertex* pv= (*it)->GetV1();
        double area= (*it)->GetArea();
        Vec3D Norm= (*it)->GetNormalVector();
        Vec3D R (pv->GetVXPos(),pv->GetVYPos(),pv->GetVZPos());
        V+=area*(R.dot(Norm,R))/3.0;
        A+=area;
    }
    m_totVolume = V;
    m_totArea = A;

}
void Analysis::UpdateEnergy(MESH &mesh)
{
    std::vector<vertex*> pV = mesh.m_pAllV;
    std::vector<links *> pL = mesh.m_pHL;
    std::vector<triangle *> pT=mesh.m_pAllT;
    Vec3D *pBox=mesh.m_pBox;

    //========== Update the Mesh geometry variables ===============
    //======== Prepare the trinagles: calculate area and normal vector
    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
        (*it)->UpdateNormal_Area(pBox);
    //===== Prepare links:  normal vector and shape operator
    for (std::vector<links *>::iterator it = pL.begin() ; it != pL.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(pBox);
    }
    //======= Prepare vertex:  area and normal vector and curvature
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
        Curvature P(*it);
    
    Energy En (m_pState->m_pinc_ForceField);
    m_totEn = En.TotalEnergy(pV,pL);

}
void Analysis::Center(MESH &mesh)
{
    std::vector<vertex*> pV = mesh.m_pAllV;
    Vec3D Box = *(mesh.m_pBox);
    Vec3D CM;

    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vec3D Pos((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        
        CM = CM+Pos*(1/double(pV.size()));
    }
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vec3D Pos((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        (*it)->NOPBCUpdatePos(Pos-CM+Box*0.5);
    }
    
    
}
void Analysis::RemovePBC(MESH &mesh)
{



    std::vector<vertex*> pV = mesh.m_pAllV;
    std::vector< L_CLUSTER > AllCluster;

    LMatrix VerMat(pV.size(),pV.size(),'I');
    for (int n=0;n<pV.size();n++)
        VerMat(n,n) = 1;
    
    

    
    for (std::vector<vertex *>::iterator it1 = pV.begin() ; it1 != pV.end(); ++it1)
    {
            int id1=(*it1)->GetVID();
            std::vector <vertex *> NV = (*it1)->GetVNeighbourVertex();
            int connect = 0;

            for (std::vector<vertex *>::iterator it3 = NV.begin() ; it3 != NV.end(); ++it3)
            {
                 connect = 0;
                int id2 = (*it3)->GetVID();

                    double dx =(*it1)->GetVXPos()-(*it3)->GetVXPos();
                    double dy =(*it1)->GetVYPos()-(*it3)->GetVYPos();
                    double dz =(*it1)->GetVZPos()-(*it3)->GetVZPos();
                    
                    
                    double dist=dx*dx+dy*dy+dz*dz;
                    //// check the distance
                    
                    if(dist<3.001)
                        connect = 1;
                if(connect == 1)
                    VerMat(id1,id2) = connect;
                }

        }


for (int n=0;n<pV.size();n++)
{
    for (int m=n+1;m<pV.size();m++)
    {
        if(VerMat(n,m)!=0)
        {
            for (int i=0;i<pV.size();i++)
            {
                if(VerMat(m,i)==0)
                    VerMat(m,i)=VerMat(n,i);
                    VerMat(n,i)=0;
                    
                    }
            break;
        }
    }
}
    
    int clid = 0 ;
    int Vsize = pV.size();
    for (int n=0;n<Vsize;n++)
    {
        int cluster = 0;
        //  std::cout<<clid<<"  "<<n<<" here 1 wrif[vm, \n";
        
        for (int m=0;m<pV.size();m++)
        {
            if(VerMat(n,m)!=0)
            {
                cluster++;
                break;
            }
        }
        
        std::vector<vertex *> TV;
        if(cluster!=0)
        {
            
            // TV.push_back(m_pAllV.at(n));
            for (int m=0;m<pV.size();m++)
            {
                if(VerMat(n,m)!=0)
                {
                    (pV.at(m))->UpdateGroup(clid);
                    //  (m_pAllV.at(m)->GetInclusion())->UpdateType("P"+f.Int_to_String(clid),clid);
                    TV.push_back(pV.at(m));
                    
                }
            }
            
            L_CLUSTER CL;
            CL.Vertex =TV;
            CL.checked = 0;
            Vec3D LX;
            CL.LX=LX;
            AllCluster.push_back(CL);
            
            clid++;
        }
    }
    std::cout<<"System contains "<<AllCluster.size()<<"  CLusters \n";
    for (std::vector<L_CLUSTER>::iterator it = AllCluster.begin() ; it != AllCluster.end(); ++it)
    {
        //  (*it).checked = 1;
        std::vector <vertex *> CV = (*it).Vertex;
        for (std::vector<vertex*>::iterator it1 = CV.begin() ; it1 != CV.end(); ++it1)
        {
            std::vector <vertex *> NV = (*it1)->GetVNeighbourVertex();
            for (std::vector<vertex*>::iterator it2 = NV.begin() ; it2 != NV.end(); ++it2)
            {
                if((*it1)->GetGroup() == (*it2)->GetGroup())
                {
                    
                }
                else
                {
                    ////====
                    int clusterid1 = (*it1)->GetGroup();
                    int clusterid2 = (*it2)->GetGroup();
                    L_CLUSTER CL2 = AllCluster.at(clusterid2);
                    bool isalready = false;
                    std::vector<int> NCluster1 = (*it).NCluster;
                    std::vector<int> NCluster2 = CL2.NCluster;
                    std::vector<vertex *> ItsV1 = (*it).ItsV;
                    std::vector<vertex *> NV1 = (*it).NV;
                    std::vector<vertex *> ItsV2 = CL2.ItsV;
                    std::vector<vertex *> NV2 = CL2.NV;
                    for (std::vector<int>::iterator itint = NCluster1.begin() ; itint != NCluster1.end(); ++itint)
                    {
                        if(clusterid2==(*itint))
                        {
                            isalready =true;
                            break;
                        }
                    }
                    if(isalready==false)
                    {
                        ItsV1.push_back(*it1);
                        NV1.push_back(*it2);
                        ItsV2.push_back(*it2);
                        NV2.push_back(*it1);
                        (*it).NV=NV1;
                        (*it).ItsV=ItsV1;
                        CL2.NV=NV2;
                        CL2.ItsV=ItsV2;
                        NCluster1.push_back(clusterid2);
                        NCluster2.push_back(clusterid1);
                        CL2.NCluster=NCluster2;
                        (*it).NCluster=NCluster1;
                    }
                }
            }
        }
    }
    
    //==
    std::vector<vertex> m_OLDV;
    
    for (std::vector<vertex*>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        m_OLDV.push_back(*(*it));
    }
    Vec3D m_OLDBox = *(mesh.m_pBox);
    //  m_Box = m_Box*10;;
    (AllCluster.at(0)).checked=1;
    
    //  Translate(m_pAllV, (m_OLDBox*5));
    
    srand(88888);
    int i=0;
    
    while (true)
    {
        int breaking = true;
        for (std::vector<L_CLUSTER>::iterator it = AllCluster.begin() ; it != AllCluster.end(); ++it)
        {
            
            //  std::cout<<i++<<"  "<<(*it).checked<<"\n";
            
            if((*it).checked==0)
            {
                breaking = false;
                break;
            }
        }
        
        if(breaking==true)
            break;
        
        
        int rng= rand() % (AllCluster.size());
        L_CLUSTER CL = AllCluster.at(rng);
        std::vector<vertex *> Vertex = CL.Vertex;
        std::vector<int> NCluster = CL.NCluster;
        std::vector<vertex *> ItsV = CL.ItsV;
        std::vector<vertex *> NV = CL.NV;
        
        bool move = false;
        int nid=0;
        int id=0;
        
        
        
        for (std::vector<int>::iterator it = NCluster.begin() ; it != NCluster.end(); ++it)
        {
            if(AllCluster.at((*it)).checked==1)
            {
                nid=(*it);
                move = true;
                break;
            }
            id++;
        }
        if(move==true && AllCluster.at(rng).checked==0)
        {
            i++;
            
            AllCluster.at(rng).checked=1;
            if(NCluster.at(id)!=nid)
                std::cout<<"error \n";
            
            int vid1=(ItsV.at(id))->GetVID();
            int vid2=(NV.at(id))->GetVID();
            Vec3D TX;
            
            
            double dx = (m_OLDV.at(vid1)).GetVXPos()-(m_OLDV.at(vid2)).GetVXPos();
            double dy = (m_OLDV.at(vid1)).GetVYPos()-(m_OLDV.at(vid2)).GetVYPos();
            double dz = (m_OLDV.at(vid1)).GetVZPos()-(m_OLDV.at(vid2)).GetVZPos();
            Vec3D DX(dx,dy,dz);
            
            
            for (int s=0;s<3;s++)
            {
                if(fabs(DX(s))<m_OLDBox(s)/2.0)
                {
                    TX(s)=0;
                }
                else
                {
                    if((DX(s))<0)
                    {
                        TX(s)=m_OLDBox(s);
                    }
                    else
                    {
                        TX(s)=-m_OLDBox(s);
                        
                    }
                }
            }
            
            
            
            Vec3D LX = (AllCluster.at((NV.at(id))->GetGroup())).LX;
            
            Vec3D LX2 = LX+TX;
            (AllCluster.at(rng)).LX=LX2;
            Translate(Vertex, LX2);
            
            
        }
        
        
    }
 
    
    


}
void Analysis::Translate(std::vector<vertex*> v, Vec3D L)
{
    for (std::vector<vertex*>::iterator it = v.begin() ; it != v.end(); ++it)
    {
        Vec3D X((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        
        Vec3D X2=X+L;
        (*it)->NOPBCUpdatePos(X2);

    }
}
void Analysis::SeparateInOut(MESH &mesh)
{

    std::vector<triangle*> T = mesh.m_pAllT;
    std::vector<vertex*> pV = mesh.m_pAllV;

    for (std::vector<vertex*>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        (*it)->UpdateGroup(0);
    }
    

    for (std::vector<triangle*>::iterator it = T.begin() ; it != T.end(); ++it)
    {
        vertex *v1= (*it)->GetV1();
        vertex *v2= (*it)->GetV2();
        vertex *v3= (*it)->GetV3();
        
        Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
        Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
        Vec3D X3(v3->GetVXPos(),v3->GetVYPos(),v3->GetVZPos());
        Vec3D X_t1=(X1+X2+X3)*(1.0/3.0);
      /* std::cout<<X_t1(0)<<" mid "<<X_t1(1)<<"\n";
              std::cout<<X1(0)<<" v1 "<<X1(1)<<"\n";
              std::cout<<X2(0)<<" v2 "<<X2(1)<<"\n";
                            std::cout<<X3(0)<<" v3 "<<X3(1)<<"\n";*/
        Vec3D N_t1=(*it)->GetNormalVector();
        bool inside=false;

        for (std::vector<triangle*>::iterator it1 = T.begin() ; it1 != T.end(); ++it1)
        {
                   // std::cout<<(*it1)->GetTriID()<<"  "<<(*it)->GetTriID()<<"   5 get to here \n";

            if((*it1)->GetTriID()!=(*it)->GetTriID())
            {
                vertex *u1= (*it1)->GetV1();
                vertex *u2= (*it1)->GetV2();
                vertex *u3= (*it1)->GetV3();
                Vec3D Y1(u1->GetVXPos(),u1->GetVYPos(),u1->GetVZPos());
                Vec3D Y2(u2->GetVXPos(),u2->GetVYPos(),u2->GetVZPos());
                Vec3D Y3(u3->GetVXPos(),u3->GetVYPos(),u3->GetVZPos());
                Vec3D X_t2=(Y1+Y2+Y3)*(1.0/3.0);
                Vec3D R=X_t2-X_t1;
                 Vec3D N_t2=(*it1)->GetNormalVector();

                double norm=R.norm();

                double cosT=R.dot(R,N_t1)*(1/norm);
                double cosTmax=norm/(4/3+norm*norm); // max angle
                
               // std::cout<<cosT<<"  "<<cosTmax<<"\n";
                if(cosT>cosTmax && N_t2.dot(N_t2,N_t1)<0)
                {
                    inside=true;
                    break;
                }
            }
        }
        if(inside==true )
        {
            v1->UpdateGroup(1);
            v2->UpdateGroup(1);
            v3->UpdateGroup(1);
        }


        
        
    }

    
}
void Analysis::SeparateInOutByVer(MESH &mesh)
{

    std::vector<triangle*> T = mesh.m_pAllT;
    std::vector<vertex*> pV = mesh.m_pAllV;

    for (std::vector<vertex*>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        (*it)->UpdateGroup(0);
    }
    

    for (std::vector<vertex*>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vec3D X_t1((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        Vec3D N_t1=(*it)->GetNormalVector();
        N_t1 = N_t1*(1.0/N_t1.norm());
        bool inside=false;

        for (std::vector<vertex*>::iterator it1 = pV.begin() ; it1 != pV.end(); ++it1)
        {
       
                   // std::cout<<(*it1)->GetTriID()<<"  "<<(*it)->GetTriID()<<"   5 get to here \n";

            if((*it1)->GetVID()!=(*it)->GetVID())
            {
       		 Vec3D X_t2((*it1)->GetVXPos(),(*it1)->GetVYPos(),(*it1)->GetVZPos());
                Vec3D R=X_t2-X_t1;
        	Vec3D N_t2=(*it1)->GetNormalVector();
                double areav2=(*it1)->GetArea();

        	N_t2 = N_t2*(1.0/N_t2.norm());
                double norm=R.norm();

                double cosT=(R.dot(R,N_t1))/norm;
                double cosTmax=norm/sqrt(areav2+norm*norm); // max angle
                

                if(cosT>cosTmax )//&& N_t2.dot(N_t2,N_t1)<-0.5 )
                {
                             //  std::cout<<cosT<<"  "<<cosTmax<<"\n";
                    inside=true;
                    break;
                }
            }
               
        }
           
        if(inside==true)
        {
            (*it)->UpdateGroup(1);
        }


     
        
    }

 
}
