#include <fstream>
#include "MESH.h"
#include "LMatrix.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
MESH class for quick access 
 */
MESH::MESH() : m_MeshCrossedPBC(false) {

    m_pEdgeV.clear();
    Vec3D b;
    m_Box = b;
    m_pBox = &m_Box;
    m_No_VectorFields_Per_V = 0;

}
MESH::~MESH() {
    
}
void MESH::RemoveFromLinkList(links* z, std::vector<links*> &vect){
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
    return;
}
void MESH::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
    return;
}
void MESH::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect){
    
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
    return;
}
void MESH::MakeALinkGhost(links* l){
    RemoveFromLinkList(l,m_pActiveL);
    RemoveFromLinkList(l,m_pHL);
    RemoveFromLinkList(l,m_pMHL);
    m_pGhostL.push_back(l);
    return;
}
void MESH::UpdateNoVFPerVertex(int number){
    
    m_No_VectorFields_Per_V = number;
    return;
}
void MESH::MakeATriangleGhost(triangle* tri){
    RemoveFromTriangleList(tri,m_pActiveT);
    m_pGhostT.push_back(tri);
    return;
}
void  MESH::CenterMesh(){
    double xcm=0;
    double ycm=0;
    double zcm=0;
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        xcm+=(*it)->GetVXPos();
        ycm+=(*it)->GetVYPos();
        zcm+=(*it)->GetVZPos();
    }
    xcm=xcm/double(m_pActiveV.size());
    ycm=ycm/double(m_pActiveV.size());
    zcm=zcm/double(m_pActiveV.size());
    
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        (*it)->UpdateVXPos((*it)->GetVXPos()-xcm+(*m_pBox)(0)/2.0);
        (*it)->UpdateVYPos((*it)->GetVYPos()-ycm+(*m_pBox)(1)/2.0);
        (*it)->UpdateVZPos((*it)->GetVZPos()-zcm+(*m_pBox)(2)/2.0);
    }
    
    return;
}

void  MESH::CenterMeshSphericalNOPBC(){
    double xcm=0;
    double ycm=0;
    double zcm=0;
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        xcm+=(*it)->GetVXPos();
        ycm+=(*it)->GetVYPos();
        zcm+=(*it)->GetVZPos();
    }
    xcm=xcm/double(m_pActiveV.size());
    ycm=ycm/double(m_pActiveV.size());
    zcm=zcm/double(m_pActiveV.size());
    
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        Vec3D Pos((*it)->GetVXPos()-xcm+(*m_pBox)(0)/2.0,(*it)->GetVYPos()-ycm+(*m_pBox)(1)/2.0,(*it)->GetVZPos()-zcm+(*m_pBox)(2)/2.0);
        (*it)->NOPBCUpdatePos(Pos);
    }
    
    return;
}

void MESH::CenterSemiFlatNOPBC()
{
    Vec3D Box = *(m_pBox);
    

    double minZ=(*m_pActiveV.begin())->GetVZPos();
    double maxZ=(*m_pActiveV.begin())->GetVZPos();


    double zpos,ypos,xpos,xavg,yavg,minX,minY;


    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){
        zpos = (*it)->GetVZPos(); // Get the Z-coordinate of the vertex
            if (zpos < minZ) {
                minZ = zpos; // Update min height
            }
            if (zpos > maxZ) {
                maxZ = zpos; // Update max height
            }

        xpos=(*it)->GetVXPos(); // Get the Z-coordinate of the vertex
        ypos=(*it)->GetVYPos(); // Get the Z-coordinate of the vertex
        xavg+=xpos;
        yavg+=ypos;
    }

    xavg=xavg/double(m_pActiveV.size());
    yavg=yavg/double(m_pActiveV.size());
    

    


    double height=(maxZ-minZ);
    double limit=minZ+height*0.125;

    minX=xavg;
    minY=yavg;

    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){

        zpos = (*it)->GetVZPos(); // Get the Z-coordinate of the vertex
            if (zpos < limit) { 
                xpos = (*it)->GetVXPos(); // Get the Z-coordinate of the vertex
                if (xpos < minX) {
                    minX = xpos; // Update min height
                }
                ypos = (*it)->GetVYPos(); // Get the Z-coordinate of the vertex
                if (ypos < minY) {
                    minY = ypos; // Update min height
                }
    }
    }

    Vec3D CM=(minX,minY,minZ);
    std::vector<vertex *> Vertex =m_pActiveV;
    //Translate(m_pActiveV,CM);
    
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        //std::cout<<"1:"<<(*it)->GetVXPos()<<"  "<<(*it)->GetVYPos()<<"  "<<(*it)->GetVZPos()<<"\n";
        Vec3D Pos((*it)->GetVXPos()-minX,(*it)->GetVYPos()-minY,(*it)->GetVZPos()-minZ);
        (*it)->NOPBCUpdatePos(Pos);
        //std::cout<<"2:"<<(*it)->GetVXPos()<<"  "<<(*it)->GetVYPos()<<"  "<<(*it)->GetVZPos()<<"\n";
    }
    
    
}


void  MESH::CenterSemiFlatMesh(){

    double minHeight=(*m_pActiveV.begin())->GetVZPos();
    double maxHeight=(*m_pActiveV.begin())->GetVZPos();;
    double zpos;


    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){

        zpos = (*it)->GetVZPos(); // Get the Z-coordinate of the vertex
            if (zpos < minHeight) {
                minHeight = zpos; // Update min height
            }
            if (zpos > maxHeight) {
                maxHeight = zpos; // Update max height
            }
    }

    double thickness=maxHeight-minHeight;
    double Lz=(*m_pBox)(2);

    if (thickness<0.9*Lz){
        return;
    }

    
    for (std::vector<vertex *>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it){

        zpos = (*it)->GetVZPos();
        
        if (zpos>0.5*Lz){
            (*it)->UpdateVZPos(zpos-0.5*Lz);
        } 
        else if (zpos<0.5*Lz){
            (*it)->UpdateVZPos(zpos+0.5*Lz);
        }

    }
    
    return;
}

struct L_CLUSTER {
    
    std::vector<vertex *> Vertex;
    std::vector<int> NCluster;
    std::vector<vertex *> ItsV;
    std::vector<vertex *> NV;
    int checked;
    Vec3D LX;
} ;

void MESH::RemovePBC()
{




    std::vector< L_CLUSTER > AllCluster;

    LMatrix VerMat(m_pActiveV.size(),m_pActiveV.size(),'I');
    for (int n=0;n<m_pActiveV.size();n++)
    {VerMat(n,n) = 1;}
    
    

    
    for (std::vector<vertex *>::iterator it1 = m_pActiveV.begin() ; it1 != m_pActiveV.end(); ++it1)
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


for (int n=0;n<m_pActiveV.size();n++)
{
    for (int m=n+1;m<m_pActiveV.size();m++)
    {
        if(VerMat(n,m)!=0)
        {
            for (int i=0;i<m_pActiveV.size();i++)
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
    int Vsize = m_pActiveV.size();
    for (int n=0;n<Vsize;n++)
    {
        int cluster = 0;
        //  std::cout<<clid<<"  "<<n<<" here 1 wrif[vm, \n";
        
        for (int m=0;m<m_pActiveV.size();m++)
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
            for (int m=0;m<m_pActiveV.size();m++)
            {
                if(VerMat(n,m)!=0)
                {
                    (m_pActiveV.at(m))->UpdateGroup(clid);
                    //  (m_pAllV.at(m)->GetInclusion())->UpdateType("P"+f.Int_to_String(clid),clid);
                    TV.push_back(m_pActiveV.at(m));
                    
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
    
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        m_OLDV.push_back(*(*it));
    }
    Vec3D m_OLDBox = *(m_pBox);
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

void MESH::Translate(std::vector<vertex*> v, Vec3D L)
{
    for (std::vector<vertex*>::iterator it = v.begin() ; it != v.end(); ++it)
    {
        Vec3D X((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        
        Vec3D X2=X+L;
        (*it)->NOPBCUpdatePos(X2);

    }
}



bool MESH::GenerateMesh(MeshBluePrint meshblueprint)
{
    m_Vertex.clear();
    m_Triangle.clear();
    m_Links.clear();
    m_Inclusion.clear();
    m_pActiveV.clear();
    m_pSurfV.clear();
    m_pEdgeV.clear();
    m_pActiveL.clear();
    m_pHL.clear();
    m_pMHL.clear();
    m_pEdgeL.clear();
    m_pActiveT.clear();
    m_pInclusion.clear();
    m_pGhostT.clear();
    m_pGhostL.clear();
    m_pGhostV.clear();
    m_Box = meshblueprint.simbox;
    if(m_Box.isbad()){
        std::cout<<"---> box from blueprint is bad \n";
        return false;
    }
    
    //m_InclusionType = meshblueprint.binctype;
    
    //for (std::vector<InclusionType>::iterator it = m_InclusionType.begin() ; it != m_InclusionType.end(); ++it)
    //    m_pInclusionType.push_back(&(*it));
    
    // Making vertices
    for (std::vector<Vertex_Map>::iterator it = (meshblueprint.bvertex).begin() ; it != (meshblueprint.bvertex).end(); ++it)
    {
            vertex v(this, it->id,it->x,it->y,it->z);
            v.UpdateBox(m_pBox);
            v.UpdateDomainID(it->domain);
            m_Vertex.push_back(v);
    }
//===== Make exclution [since June, 2023]
       
       std::vector<int> excluded_ver = meshblueprint.excluded_id; // this vertices should be excluded
       for (std::vector<int>::iterator it = excluded_ver.begin() ; it != excluded_ver.end(); ++it)
       {
           ((meshblueprint.bvertex).at((*it))).include = false;
       }
       int t=0;
       for (std::vector<vertex>::iterator it = m_Vertex.begin() ; it != m_Vertex.end(); ++it)
       {
               if(((meshblueprint.bvertex).at((t))).include == true)
                   m_pActiveV.push_back(&(*it));
           t++;
       }
    //== since the vertex number remain constant, it is better to change the id of active vertices
    t=0;
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        (*it)->UpdateVID(t);
        t++;
    }
    // Making triangles
    t=0;
    int temid = 0;
    for (std::vector<Triangle_Map>::iterator it = (meshblueprint.btriangle).begin() ; it != (meshblueprint.btriangle).end(); ++it)
    {
        bool pr=true;
        
        int vid1 = ((meshblueprint.btriangle).at(t)).v1;
        int vid2 = ((meshblueprint.btriangle).at(t)).v2;
        int vid3 = ((meshblueprint.btriangle).at(t)).v3;

        if(((meshblueprint.bvertex).at(vid1)).include == true)
        if(((meshblueprint.bvertex).at(vid2)).include == true)
        if(((meshblueprint.bvertex).at(vid3)).include == true)
        {
        triangle T(temid,&(m_Vertex.at(it->v1)),&(m_Vertex.at(it->v2)),&(m_Vertex.at(it->v3)));
        m_Triangle.push_back(T);
            temid++;
        }
        t++;
    }
    //make inclusions
    for (std::vector<Inclusion_Map>::iterator it = (meshblueprint.binclusion).begin() ; it != (meshblueprint.binclusion).end(); ++it)
    {
        
        if(m_Vertex.size()<it->vid+1 && it->tid > m_InclusionType.size()) {
            
            std::cout<<"----> Error: Inclusion vertex id or type id is out of range "<<std::endl;
            exit(0);
        }
        inclusion Tinc(it->id, &(m_InclusionType[it->tid]));
        Tinc.Updatevertex(&(m_Vertex.at(it->vid)));
        Vec3D D(it->x,it->y,0);
        Tinc.UpdateLocalDirection(D);
        m_Inclusion.push_back(Tinc);
        (m_Vertex.at(it->vid)).UpdateOwnInclusion(true);
    }
    for (std::vector<inclusion>::iterator it = m_Inclusion.begin() ; it != m_Inclusion.end(); ++it)
        m_pInclusion.push_back(&(*it));


    
    t=0;
    for (std::vector<triangle>::iterator it = m_Triangle.begin() ; it != m_Triangle.end(); ++it)
        m_pActiveT.push_back(&(*it));

// end Make exclution
    
    
    int li=-1;
    
    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
    {
        ((*it)->GetV1())->AddtoTraingleList((*it));
        ((*it)->GetV1())->AddtoNeighbourVertex(((*it)->GetV2()));
        ((*it)->GetV2())->AddtoTraingleList((*it));
        ((*it)->GetV2())->AddtoNeighbourVertex(((*it)->GetV3()));
        ((*it)->GetV3())->AddtoTraingleList((*it));
        ((*it)->GetV3())->AddtoNeighbourVertex(((*it)->GetV1()));
        
        /// create links
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        
        links l1(id1,(*it)->GetV1(),(*it)->GetV2(),(*it));
        l1.UpdateV3((*it)->GetV3());
        
        links l2(id2,(*it)->GetV2(),(*it)->GetV3(),(*it));
        l2.UpdateV3((*it)->GetV1());
        
        links l3(id3,(*it)->GetV3(),(*it)->GetV1(),(*it));
        l3.UpdateV3((*it)->GetV2());
        m_Links.push_back(l1);
        m_Links.push_back(l2);
        m_Links.push_back(l3);
        
    }
    li=-1;
    for (std::vector<triangle*>::iterator it = m_pActiveT.begin() ; it != m_pActiveT.end(); ++it)
    {
        li++;
        int id1=li;
        li++;
        int id2=li;
        li++;
        int id3=li;
        links * l1=&(m_Links.at(id1));
        links * l2=&(m_Links.at(id2));
        links * l3=&(m_Links.at(id3));
        l1->UpdateNeighborLink1(l2);
        l1->UpdateNeighborLink2(l3);
        l2->UpdateNeighborLink1(l3);
        l2->UpdateNeighborLink2(l1);
        l3->UpdateNeighborLink1(l1);
        l3->UpdateNeighborLink2(l2);
        
        
        ((*it)->GetV1())->AddtoLinkList(l1);
        ((*it)->GetV2())->AddtoLinkList(l2);
        ((*it)->GetV3())->AddtoLinkList(l3);
        
    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        bool foundM=false;
        if((it)->GetMirrorFlag()==true)
        {
            m_pMHL.push_back(it->GetMirrorLink());
            m_pHL.push_back(&(*it));
            foundM = true;
        }
        else
        {
            vertex *v1=it->GetV1();
            vertex *v2=it->GetV2();
            
            std::vector<links*>  lList = v2->GetVLinkList();
            for (std::vector<links*>::iterator it2 = lList.begin() ; it2 != lList.end(); ++it2)
            {
                if(((*it2)->GetV2())->GetVID()==v1->GetVID())
                {
                    it->UpdateMirrorLink((*it2));
                    (*it2)->UpdateMirrorLink(&(*it));
                    it->UpdateMirrorFlag(true);
                    (*it2)->UpdateMirrorFlag(true);
                    foundM = true;
                    break;
                }
            }
        }
        if(foundM == false)
        {
            m_pEdgeL.push_back(&(*it));
            it->m_LinkType = 1;
        }
        
    }
    int edgelink=0;
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
        edgelink++;

    if(edgelink!=0)
    {
        std::cout<<"----> Note: "<<edgelink<<" links at the edge \n";
        std::cout<<"----> Note (Warning): the system is not closed! \n";

    }
    for (std::vector<links>::iterator it = m_Links.begin() ; it != m_Links.end(); ++it)
    {
        
        m_pActiveL.push_back(&(*it));
    }
    //==== Getting the edge vertex from link
    for (std::vector<links*>::iterator it = m_pEdgeL.begin() ; it != m_pEdgeL.end(); ++it)
    {
        m_pEdgeV.push_back((*it)->GetV1());
        ((*it)->GetV1())->m_pEdgeLink = *it;
        ((*it)->GetV2())->m_pPrecedingEdgeLink = *it;
        ((*it)->GetV2())->AddtoNeighbourVertex((*it)->GetV1());
        ((*it)->GetV1())->m_VertexType = 1;
    }
    for (std::vector<vertex*>::iterator it = m_pActiveV.begin() ; it != m_pActiveV.end(); ++it)
    {
        bool is_an_edge = false;
        for (std::vector<vertex*>::iterator it2 = m_pEdgeV.begin() ; it2 != m_pEdgeV.end(); ++it2) {
            if( *it2 == *it)
                is_an_edge = true;
        }
        if(!is_an_edge)
        m_pSurfV.push_back(*it);
    }
    
    //
    
    
    for (std::vector<inclusion*>::iterator it = m_pInclusion.begin() ; it != m_pInclusion.end(); ++it) {
        
        ((*it)->Getvertex())->UpdateInclusion((*it));
    }
    // =======
    // ==== info of the mesh
    
    //====== create ghost
    int lid = 2*(m_pMHL.size())+m_pEdgeL.size();
    int tid = m_pActiveT.size() ;
    int vid = m_pActiveV.size() ;

    for (int i = 0; i < m_pActiveV.size(); ++i) {
        m_GhostL.push_back(links(lid++));
    }
    for (int i = 0; i < m_pActiveV.size(); ++i) {
        m_GhostT.push_back(triangle(tid++));
    }
    for (int i = 0; i < m_pActiveV.size(); ++i) {
        m_GhostV.push_back(vertex(this, vid++));
    }

    for (std::vector<links>::iterator it = m_GhostL.begin() ; it != m_GhostL.end(); ++it)
        m_pGhostL.push_back(&(*it));
    
    for (std::vector<triangle>::iterator it = m_GhostT.begin() ; it != m_GhostT.end(); ++it)
        m_pGhostT.push_back(&(*it));
    
    for (std::vector<vertex>::iterator it = m_GhostV.begin() ; it != m_GhostV.end(); ++it)
        m_pGhostV.push_back(&(*it));

    
    //--- setting up the vector fields
    m_No_VectorFields_Per_V = meshblueprint.number_vector_field;

    // Ensure the size of meshblueprint.bvectorfields matches the number of active vertices
    if (m_No_VectorFields_Per_V !=0 && meshblueprint.bvectorfields.size() != m_pActiveV.size()) {
        std::cerr << "Error: Mismatch between vector fields data and active vertices count.\n";
        std::cerr << "---- active vertices: "<<m_pActiveV.size()<<"  vector field data "<<meshblueprint.bvectorfields.size()<<"\n";
        return false;
    }

    if(m_No_VectorFields_Per_V != 0 ) {
        // Get an iterator to the beginning of meshblueprint.bvectorfields
        std::cout<<"---> Note, each vertex has "<< m_No_VectorFields_Per_V <<" vector fields \n";
        std::vector<VectorField_Map>::iterator data_it = meshblueprint.bvectorfields.begin();
        // Iterate over active vertices and initialize them with corresponding vector field data
        for (std::vector<vertex*>::iterator it = m_pActiveV.begin(); it != m_pActiveV.end(); ++it, ++data_it) {
            std::string data_line = data_it->data_line;
            (*it)->Initialize(m_No_VectorFields_Per_V, data_line, this);
        }
    }//     if(m_No_VectorFields_Per_V !=0) {
    else{
    }
    
    if(m_No_VectorFields_Per_V != 0 )
    {
        for (std::vector<links *>::const_iterator it = m_pActiveL.begin() ; it != m_pActiveL.end(); ++it) {
            (*it)->InitializeVFIntEnergy(m_No_VectorFields_Per_V);
        }
        for (std::vector<links *>::const_iterator it = m_pGhostL.begin() ; it != m_pGhostL.end(); ++it) {
            (*it)->InitializeVFIntEnergy(m_No_VectorFields_Per_V);
        }
    }
    //WritevtuFiles VTU(pState);
    //std::string file="ini_Mesh.vtu";
    //VTU.Writevtu(m_pAllV,m_pAllT,m_pAllLinks,file);
    
    return true;
}
//===========================================================
// Note, the converted blue print will not have the exclusions
//
MeshBluePrint MESH::Convert_Mesh_2_BluePrint(MESH *mesh) {
    
    
    MeshBluePrint BluePrint;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector <InclusionType> binctype;  // a vector containing all inclsuion type and a default one
    std::vector <VectorField_Map> bvectorfields;  // a vector containing all inclsuion type and a default one

    Vec3D simbox;
    
    // vertex member of the blue print
    std::vector<vertex*> pV = mesh->m_pActiveV;
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        Vertex_Map tvm;
        tvm.x = (*it)->GetVXPos();
        tvm.y = (*it)->GetVYPos();
        tvm.z = (*it)->GetVZPos();
        tvm.id = (*it)->GetVID();
        tvm.domain = (*it)->GetGroup();
        bvertex.push_back(tvm);
        
        //--- vector fields
        VectorField_Map tem_vf;
        tem_vf.data_line = (*it)->GetVectorFieldsStream();
        bvectorfields.push_back(tem_vf);
    }
    // triangle map member of the blue print
    std::vector<triangle*> pT = mesh->m_pActiveT;
    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
    {
        Triangle_Map ttm;
        ttm.v1 = ((*it)->GetV1())->GetVID();
        ttm.v2 = ((*it)->GetV2())->GetVID();
        ttm.v3 = ((*it)->GetV3())->GetVID();
        ttm.id = (*it)->GetTriID();
        btriangle.push_back(ttm);
    }
    
    // inclusion map member of the blue print
    std::vector<inclusion*> pInc = mesh->m_pInclusion;
    for (std::vector<inclusion *>::iterator it = pInc.begin() ; it != pInc.end(); ++it)
    {
        Inclusion_Map tim;
        tim.x = ((*it)->GetLDirection())(0);
        tim.y = ((*it)->GetLDirection())(1);
        tim.vid = (*it)->Getvertex()->GetVID();
        tim.tid = (*it)->GetInclusionType()->ITid;
        tim.id = (*it)->GetID();
        binclusion.push_back(tim);

    }
    //--- note vector fields are written in the vertex section
    
    // Add other map into the mesh map
    BluePrint.bvertex = bvertex;
    BluePrint.btriangle = btriangle;
    BluePrint.binclusion = binclusion;
    BluePrint.bvectorfields = bvectorfields;
    BluePrint.number_vector_field = mesh->GetNoVFPerVertex();
    BluePrint.simbox = *(mesh->m_pBox);
    
    return BluePrint;
}
bool MESH::UpdateGroupFromIndexFile(std::string &filename){
    //std::map<std::string, std::vector<vertex*>> GetGroups()  const  {return m_Groups;}

    if (Nfunction::SubstringFromRight(filename, '.') != "inx") {
        filename += ".inx";
    }
    
    std::ifstream indexfile(filename);
    if (!indexfile) {
        std::cerr << "---> note: no index file has been provided "  << std::endl;
        return false;
    }

    std::string name;
    int NAtom;
    int gid = 1;
    while (indexfile >> name >> NAtom) {
        int vid;
        std::vector<vertex*> gmember;
        for (int i = 0; i < NAtom; ++i) {
            if (!(indexfile >> vid)) {
                std::cerr << "---> error: Failed to read ID from index file" << std::endl;
                return false;
            }
            if (vid < m_pActiveV.size()) {
                m_pActiveV[vid]->UpdateGroupName(name);
                m_pActiveV[vid]->UpdateGroup(gid);
                gmember.push_back(m_pActiveV[vid]);
            } else {
                std::cerr << "---> error: ID " << vid << " exceeds active vector size" << std::endl;
                return false;
            }
        }
        m_Groups[name] = gmember;
        gid++;
    }
    indexfile.close();
    return true;
}
double MESH::SquareDistanceBetweenTwoVertices(vertex *p_v1, vertex* p_v2, Vec3D Box){
        
    double dx = p_v1->GetVXPos() - p_v2->GetVXPos();
    double dy = p_v1->GetVYPos() - p_v2->GetVYPos();
    double dz = p_v1->GetVZPos() - p_v2->GetVZPos();

    double boxHalfX = Box(0) / 2.0;
    double boxHalfY = Box(1) / 2.0;
    double boxHalfZ = Box(2) / 2.0;

    // Adjust coordinates if outside the periodic boundary
    if (fabs(dx) > boxHalfX) {
        dx = (dx < 0) ? Box(0) + dx : dx - Box(0);
    }
    if (fabs(dy) > boxHalfY) {
        dy = (dy < 0) ? Box(1) + dy : dy - Box(1);
    }
    if (fabs(dz) > boxHalfZ) {
        dz = (dz < 0) ? Box(2) + dz : dz - Box(2);
    }

    // Compute and return squared distance
    return dx * dx + dy * dy + dz * dz;
    
}


