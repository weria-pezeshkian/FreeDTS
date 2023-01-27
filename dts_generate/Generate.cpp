

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Generate.h"
Generate::Generate()
{
}
Generate::Generate(std::vector <std::string> argument)
{

    m_Argument = argument;
    m_Healthy =true;
    m_Seed =36723;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_OutputFilename = "out.q";
    m_Box(0) = 10;
    m_Box(1) = 10;
    m_Box(2) = 10;
    m_genus = 0;
    m_Type = "Flat";
    m_N = 5;
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);

    ExploreArguments();     // read the input data


    // After reading the argument we generate a TS file based on the provided values
    if (m_Type == "flat" || m_Type == "Flat" || m_Type == "FLAT" )
    {
        FlatBilayer();
    }
    else if (m_Type == "Tetrahedron" || m_Type == "tetrahedron"  )
    {
        Tetrahedron();
    }
    else if (m_Type == "channel" || m_Type == "Channel"  )
    {
        Cylinder();
    }
    else if (m_Type == "high_gen" || m_Type == "High_gen"  )
    {
        HighTopologyStructure();
    }
    else
    {
        std::cout<<" Shape type "<<m_Type<<" is not recognized "<<std::endl;
    }
    
}
Generate::~Generate()
{
    
}
//======== high genus
void Generate::HighTopologyStructure()
{
    srand (40072);


///=======
//== Read Gen-morphology variables
//===========
    int TopDegree = m_genus;
    int NoVertex = m_N;
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    double DR = 0.001;
    double dx = 0;
               // cleans the log and error files
    
    if(m_genus>m_N*m_N)
    {
        std::cout<<" the topology degree is high, you may use higher number of verteces \n";
    }
    m_noUpperCreatedHole = 0;
    m_noLowerCreatedHole = 0;
//============================
// Gen vertex
//==============================
    double zm = m_Box(2)/2;
    double xm= 1;
    double ym= 1;
    double l=1.01;
    int id=0;
    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
         
            
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm+l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
        
    }
    
    

    for (int i=0;i<m_N;i++)
    {
        for (int j=0;j<m_N;j++)
        {
            
            double x = i*l+xm;
            double y = j*l+ym;
            double z = zm-l/2;
            Vertex_Map v;
            v.id = id;
            dx = double(rand()%1000)/1000.0;
            v.x = x+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.y = y+dx*DR;
            dx = double(rand()%1000)/1000.0;
            v.z = z+dx*DR;
            v.domain = 0;
            allV.push_back(v);
            id++;
        }
        
    }
    

//============================
// Connect triangles
//==============================
   id=0;
    bool makehole=false;
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            

            makehole=MakeHole(1,i,j);
            
            if(makehole==true)
            {
                std::vector<Triangle_Map> temT = MakeTrianglesAroundHole(id, i, j);
                for (std::vector<Triangle_Map>::iterator it = temT.begin() ; it != temT.end(); ++it)
                    allT.push_back(*it);

                id=id+8;
            }
            else
            {
                int v1id=idfromij(1,i,j);
                int v2id=idfromij(1,i,j+1);
                int v3id=idfromij(1,i+1,j+1);
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                
                v1id=idfromij(1,i,j);
                v2id=idfromij(1,i+1,j+1);
                v3id=idfromij(1,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }


            
        }
        
    }
    
    // lower layer
    for (int i=0;i<m_N-1;i++)
    {
        for (int j=0;j<m_N-1;j++)
        {
            int v1id=idfromij(2,i,j);
            int v3id=idfromij(2,i,j+1);
            int v2id=idfromij(2,i+1,j+1);

            makehole=MakeHole(2,i,j);

            
            if(makehole==true)
            {

            }
            else
            {
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                v1id=idfromij(2,i,j);
                v3id=idfromij(2,i+1,j+1);
                v2id=idfromij(2,i+1,j);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }
            

            
        }
        
    }
    

    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,0,i);
        int v2id=idfromij(2,0,i);
        int v3id=idfromij(2,0,i+1);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
         v1id=idfromij(1,0,i);
         v2id=idfromij(2,0,i+1);
         v3id=idfromij(1,0,i+1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;

        
    }
    
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,m_N-1,i+1);
        int v2id=idfromij(2,m_N-1,i+1);
        int v3id=idfromij(2,m_N-1,i);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,m_N-1,i+1);
        v2id=idfromij(2,m_N-1,i);
        v3id=idfromij(1,m_N-1,i);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,i+1,0);
        int v2id=idfromij(2,i+1,0);
        int v3id=idfromij(2,i,0);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,i+1,0);
        v2id=idfromij(2,i,0);
        v3id=idfromij(1,i,0);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    
    
    for (int i=0;i<m_N-1;i++)
    {
        
        int v1id=idfromij(1,i,m_N-1);
        int v2id=idfromij(2,i,m_N-1);
        int v3id=idfromij(2,i+1,m_N-1);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,i,m_N-1);
        v2id=idfromij(2,i+1,m_N-1);
        v3id=idfromij(1,i+1,m_N-1);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }
    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
int Generate::idfromij(int s, int i, int j)
{
    
    if(i>=m_N || m_N<=j)
    {
        std::cout<<"error, this should happen \n";
    }
    int id=0;
    if(s==1)
    id=m_N*j+i;
    else if(s==2)
    id=m_N*m_N+m_N*j+i;
    
    
    
    return id;
}
std::vector<Triangle_Map>  Generate::MakeTrianglesAroundHole(int id, int i, int j)
{
    
    std::vector<Triangle_Map> allT;

    int v1id=idfromij(1,i,j);
    int v2id=idfromij(1,i,j+1);
    int v3id=idfromij(2,i,j+1);
    Triangle_Map T1;
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;


    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i,j+1);
    v3id=idfromij(2,i,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    
    v1id=idfromij(1,i,j+1);
    v2id=idfromij(2,i+1,j+1);
    v3id=idfromij(2,i,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j+1);
    v2id=idfromij(1,i+1,j+1);
    v3id=idfromij(2,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i,j);
    v3id=idfromij(2,i+1,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j);
    v2id=idfromij(2,i+1,j);
    v3id=idfromij(1,i+1,j);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j);
    v2id=idfromij(2,i+1,j+1);
    v3id=idfromij(1,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j);
    v2id=idfromij(2,i+1,j);
    v3id=idfromij(2,i+1,j+1);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);

    return allT;
}

bool Generate::MakeHole(int s,int i, int j)
{
    bool is=false;
    
    if(m_genus!=0)
    {
    
    int space=(m_N-2)/(int(sqrt(m_genus))+1);
    

        if(s==1 && m_noUpperCreatedHole<m_genus)
        {
            
            if((i+1)%space==0 && (j+1)%space==0 && j<m_N-1 && i<m_N-1)
                is=true;
            
        }
        else if(s==2 && m_noLowerCreatedHole<m_genus && j<m_N-1 && i<m_N-1)
        {
            if((i+1)%space==0 && (j+1)%space==0)
                is=true;
        }
    }



    
    
    
    if(is==true && s==1)
    m_noUpperCreatedHole++;
    else if(is==true && s==2)
    m_noLowerCreatedHole++;
    
    return is;
}
//======= end high genus
void Generate::Cylinder()
{
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    
    double l=1.2;
    int Nz = int(m_Box(2)/l);
    m_Box(2) = double(Nz)*l;
    double lx=l*0.87;
    
    int Nx= m_N;
    if (Nx%2!=0)
        Nx=Nx+1;
    double Width = double(Nx)*lx;
    
    double xm=m_Box(1)/2.0;
    
    int id=0;
    
    for (int i=0;i<Nz;i++)
    {
        for (int j=0;j<Nx;j++)
        {
            
            double y=(double(j)+0.5)*lx-Width/2+xm;
            double x=xm+Width/2+0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
       }
        for (int j=0;j<Nx;j++)
        {
            
            double x=-(double(j)+0.5)*lx+Width/2+xm;
            double y=xm+Width/2+0.2;
            double z=l*(double(i)+double((j+Nx)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
        for (int j=0;j<Nx;j++)
        {
            
            double y=-(double(j)+0.5)*lx+Width/2+xm;
            double x=-Width/2+xm-0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
        for (int j=0;j<Nx;j++)
        {
            
            double x=(double(j)+0.5)*lx-Width/2+xm;
            double y=-Width/2+xm-0.2;
            double z=l*(double(i)+double((j)%2)*0.5);
            
            Vertex_Map v;
            v.id = id;
            v.x = x;
            v.y = y;
            v.z = z;
            v.domain = 0;
            allV.push_back(v);
            
            id++;
            
        }
    }
 
// Making triangles
    int t=0;
    
    for (int i=0;i<Nz;i++)
    {
        for(int j=0;j<4*Nx;j++)
        {
            int M = 4*Nx;
            if(j%2==0)
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j);
                allT.push_back(T1);
                
                t++;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i+1,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j+1);
                allT.push_back(T1);
                t++;
                
            }
            else
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j+1);
                allT.push_back(T1);
                t++;
                
                T1.id = t;
                T1.v1 = CylinderIndex(M,Nz,i,j);
                T1.v2 = CylinderIndex(M,Nz,i+1,j+1);
                T1.v3 = CylinderIndex(M,Nz,i+1,j);
                allT.push_back(T1);
                t++;
                
                
            }
        }
    }
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
void Generate::Tetrahedron()
{
    
    int N=m_N;
    int Nv=2*(N*N+1);
    
    
    double x=0;
    double y=0;
    double z=0;
    int id=0;
    double b=1.2;
    double dl=N*b;
    double theta=asin(1.0/sqrt(3.0));
    double H=dl*cos(theta);
    int M;
    
    
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
  
    double lx=m_Box(0)/2;
    double ly=m_Box(1)/2;
    double lz=m_Box(2)/2;
    
    
    //novertex=3*(N*N)+2
    for (int i=0;i<N+1;i++)
    {
        if(i<1)
            M=0;
        else if(i>=1)
            M=i-1;
        for (int j=0;j<M+1;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*i*sin(theta)/2.0;
            y=-dl/2.0+j*b+(N-i)*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
        }
    }
    for (int i=0;i<N+1;i++)
    {
        for (int j=0;j<i;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*(i-j)*sqrt(3.0)/2.0-b*i*sin(theta);
            y=b*i/2.0-j*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }
    
    for (int i=0;i<N+1;i++)
    {
        for (int j=0;j<i;j++)
        {
            
            z=H-b*i*cos(theta);
            x=b*j*sqrt(3.0)/2.0-b*i*sin(theta);
            y=-j*b/2.0;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }
    
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i;j++)
        {
            
            z=0;
            double y1=dl*sqrt(3)/3-i*b*sqrt(3)/2;
            double x1=-i*b/2+j*b;
            x=x1;
            y=y1;
            double C=sqrt(3)/2;
            double S=-0.5;
            x=C*x1-S*y1;
            y=S*x1+C*y1;
            Vertex_Map v;
            v.id = id;
            v.x = x+lx;
            v.y = y+ly;
            v.z = z+lz;
            v.domain = 0;
            allV.push_back(v);
            id++;
            
        }
    }

    int Tid=0;
    
    //==================== adding three triangles associated with upper vertex
    int v1=0;
    int v2=1;
    int v3=N*(N+1)/2+1;
    
    Triangle_Map T1;
    T1.id = Tid;
    T1.v1 = v1;
    T1.v2 = v2;
    T1.v3 = v3;
    allT.push_back(T1);
    Tid++;
    
    v2=N*(N+1)/2+1;
    v3=N*(N+1)+1;
    Triangle_Map T2;
    T2.id = Tid;
    T2.v1 = v1;
    T2.v2 = v2;
    T2.v3 = v3;
    allT.push_back(T2);
    Tid++;
    
    v2=N*(N+1)+1;
    v3=1;

    Triangle_Map T3;
    T3.id = Tid;
    T3.v1 = v1;
    T3.v2 = v2;
    T3.v3 = v3;
    allT.push_back(T3);
    Tid++;
    //========================

    id=1;
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(0,i,j);
            v2=findid(0,i+1,j);
            v3=findid(0,i+1,j+1);

            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            
            Tid++;
            
            v1=findid(0,i,j);
            v2=findid(0,i+1,j+1);
            v3=findid(0,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            if(j==1)
            {
                
                v1=findid(0,i,j);
                v2=findid(0,i+1,j-1);
                v3=findid(0,i+1,j);
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    
    
  
    
    id=1;
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(1,i,j);
            v2=findid(1,i+1,j);
            v3=findid(1,i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findid(1,i,j);
            v2=findid(1,i+1,j+1);
            v3=findid(1,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            if(j==1)
            {
                
                v1=findid(1,i,j);
                v2=findid(1,i+1,j-1);
                v3=findid(1,i+1,j);
                Triangle_Map T;
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    
 
    for (int i=1;i<N;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findid(2,i,j);
            v2=findid(2,i+1,j);
            v3=findid(2,i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findid(2,i,j);
            v2=findid(2,i+1,j+1);
            v3=findid(2,i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            
            if(j==1)
            {
                
                v1=findid(2,i,j);
                v2=findid(2,i+1,j-1);
                v3=findid(2,i+1,j);
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }
    //  three extra traingales
    
    {
    v2=findid(0,N,1);
    v1=findid(0,N,2);
    v3=findid(2,N,N);
    Triangle_Map T;
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    v2=findid(1,N,1);
    v1=findid(1,N,2);
    v3=findid(0,N,N);
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    v2=findid(2,N,1);
    v1=findid(2,N,2);
    v3=findid(1,N,N);
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    }
    
    // last face a hard one
    for (int i=1;i<N-1;i++)
    {
        for (int j=1;j<i+1;j++)
        {
            v1=findidface4(i,j);
            v3=findidface4(i+1,j);
            v2=findidface4(i+1,j+1);
            Triangle_Map T;
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            v1=findidface4(i,j);
            v3=findidface4(i+1,j+1);
            v2=findidface4(i,j+1);
            T.id = Tid;
            T.v1 = v1;
            T.v2 = v2;
            T.v3 = v3;
            allT.push_back(T);
            Tid++;
            
            
            if(j==i)
            {
                
                v1=findidface4(i-1,j+1);
                v2=findidface4(i,j+1);
                v3=findidface4(i,j);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
            if(j==1)
            {
                
                v1=findidface4(i-1,j);
                v2=findidface4(i,j);
                v3=findidface4(i-1,j-1);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
                
                
                v1=findidface4(i-1,j-1);
                v2=findidface4(i,j);
                v3=findidface4(i,j-1);
                
                T.id = Tid;
                T.v1 = v1;
                T.v2 = v2;
                T.v3 = v3;
                allT.push_back(T);
                Tid++;
            }
        }
    }

    
    
    v1=findidface4(N-2,1);
    v2=findid(2,N,2);
    v3=findid(1,N,N);
    
    
    Triangle_Map T;
    T.id = Tid;
    T.v1 = v1;
    T.v2 = v2;
    T.v3 = v3;
    allT.push_back(T);
    Tid++;
    
    
    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = m_Box;
    
    std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    
    if(ext==TSExt)
    {
        WriteQFile(m_OutputFilename , BluePrint);
    }
    else if(ext==TSIExt)
    {
        WriteTSI(m_OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }
    
}
void Generate::FlatBilayer() // making flat bilayers
{
    
    srand (800);
    double lx=sqrt(2)/sqrt(2);
    double ly=(sqrt(6.0)/2.0)/sqrt(2);
    lx=lx*1.2;
    ly=ly*1.2;
    
    
    int N=int (m_Box(1)/ly);
    int M=int (m_Box(0)/lx);
    
    lx=m_Box(0)/double(M);
    ly=m_Box(1)/double(N);

    
    int t=0;
    m_Box(0)=lx*double(M);
    m_Box(1)=ly*double(N);
    
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;

// making the vertices
    for (int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        {
            double s1=0.4*(double(rand()%2000000)/2000000-0.5);
            double s2=0.4*(double(rand()%2000000)/2000000-0.5);
            double x=((i)%2*0.5+double(j))*lx;
            double y=(double(i))*ly;
            double z=m_Box(2)/2+0.02*(2+sin(y)+cos(2*x));
            
            Vertex_Map v;
            v.x = x;
            v.y = y;
            v.z = z;
            v.id = t;
            v.domain = 0;
            allV.push_back(v);
            t++;
        }
    }

    t=0;
    for (int j=0;j<N;j++)
    {
        for(int i=0;i<M;i++)
        {
            if(j%2==0)
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = findindex(M,N,i,j);
                T1.v2 = findindex(M,N,i,j+1);
                T1.v3 = findindex(M,N,i-1,j+1);
                allT.push_back(T1);
                t++;
                Triangle_Map T2;
                T2.id = t;
                T2.v1 = findindex(M,N,i,j);
                T2.v2 = findindex(M,N,i+1,j);
                T2.v3 = findindex(M,N,i,j+1);
                allT.push_back(T2);
                t++;
                
            }
            else
            {
                Triangle_Map T1;
                T1.id = t;
                T1.v1 = findindex(M,N,i,j);
                T1.v2 = findindex(M,N,i+1,j+1);
                T1.v3 = findindex(M,N,i,j+1);
                allT.push_back(T1);
                t++;
                Triangle_Map T2;
                T2.id = t;
                T2.v1 = findindex(M,N,i,j);
                T2.v2 = findindex(M,N,i+1,j);
                T2.v3 = findindex(M,N,i+1,j+1);
                allT.push_back(T2);
                t++;
            }
        }
    }
        
        MeshBluePrint BluePrint;
        BluePrint.bvertex = allV;
        BluePrint.btriangle = allT;
        BluePrint.binclusion = allI;
        BluePrint.simbox = m_Box;
        
        std::string ext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
        
        if(ext==TSExt)
        {
            WriteQFile(m_OutputFilename , BluePrint);
        }
        else if(ext==TSIExt)
        {
            WriteTSI(m_OutputFilename , BluePrint);
        }
        else
        {
            std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
        }
        
}
void Generate::WriteQFile(std::string filename , MeshBluePrint blueprint)
{
    int pres=10;
    std::ofstream output;
    output.open(filename.c_str());
    
    
    output<<std::fixed;
    output<<std::setprecision( pres )<<(blueprint.simbox)(0)<<"   "<<(blueprint.simbox)(1)<<"   "<<(blueprint.simbox)(2)<<"   \n";
    output<<(blueprint.bvertex).size()<<"\n";
    
    
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        output<<std::setprecision( pres )<<it->id<<"  "<<it->x<<"  "<<it->y<<"  "<<it->z<<"  "<<it->domain<<std::endl;

    output<< (blueprint.btriangle).size()<<"\n";
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        output<<it->id<<"  "<<it->v1<<"   "<<it->v2<<"  "<<it->v3<<"  0  "<<std::endl;
    
    output.close();
}
void Generate::WriteTSI(std::string filename , MeshBluePrint blueprint)
{

    FILE * output;
    output = fopen(filename.c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
    //------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(blueprint.simbox)(0),(blueprint.simbox)(1),(blueprint.simbox)(2));
    
    const char* ver="vertex";
    int size=(blueprint.bvertex).size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->x,it->y,it->z);
  
    const char* tri="triangle";
    size = (blueprint.btriangle).size();
    fprintf(output,"%s%20d\n",tri,size);
   for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",it->id,it->v1,it->v2,it->v3);
    
 
    const char* inc="inclusion";
    size = (blueprint.binclusion).size();
    fprintf(output,"%s%20d\n",inc,size);
     format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
     for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->tid,it->vid,it->x,it->y);

    fclose(output);
}
void Generate::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-o")
        {
            m_OutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-seed")
        {
            m_Seed = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-g")
        {
            m_genus = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-N")
        {
            m_N = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-box")
        {
            m_Box(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Box(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Box(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;
            i++;

        }
        else if(Arg1=="-type")
        {
            m_Type = m_Argument.at(i+1);
        }
        else
        {
            std::cout << "---> Error: "<<Arg1;
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }

}
void Generate::HelpMessage()
{
    std::cout<<"--------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"---------------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting  -------------------"<<"\n";
    std::cout<<" ./GEN -box 30 30 30 -type channel -N 10  -o out.tsi"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -Box     3*double       10 10 10          box sides "<<"\n";
    std::cout<<"  -o         string       out.q             out put file name, could be both tsi or q file formats by giving the extension "<<"\n";
    std::cout<<"  -type      string       flat              ts shape (flat/tetrahedron/channel/high_gen)  "<<"\n";
    std::cout<<"  -N         int          5                 vertex per side for channel and tetrahedron shape "<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
}
int Generate::findindex(int M,int N,int i,int j)
{
    int q=(i+M)%M;
    int p=j%N;
    int s=q+M*p;
    return s;
}
int Generate::CylinderIndex(int M,int N,int i,int j)
{
    int q=(j+M)%M;
    int p=(i+N)%N;
    int s=q+M*p;
    return s;
}
int Generate::findid(int faceno,int i, int j)
{
    int id=0;
    
    
    if(j>0 && i>=j && i<=m_N+1)
    {
        id=i*(i-1)/2+j+faceno*m_N*(m_N+1)/2;
    }
    else if(j>i && i<=m_N+1)
    {
        if(faceno==0)
        {
            id=i*(i-1)/2+1+m_N*(m_N+1)/2;
        }
        if(faceno==1)
        {
            id=i*(i-1)/2+1+m_N*(m_N+1);
        }
        if(faceno==2)
        {
            id=i*(i-1)/2+1;
        }
        
    }
    else if(j==0)
    {
        if(faceno==0)
        {
            id=i*(i-1)/2+i+m_N*(m_N+1);
        }
        if(faceno==1)
        {
            id=i*(i-1)/2+i;
        }
        if(faceno==2)
        {
            id=i*(i-1)/2+i+m_N*(m_N+1)/2;
        }
        
    }
    
    return id;
    
}
int Generate::findidface4(int i, int j)
{
    int id=0;
    
    if(j>0 && i>=j && i<=m_N-2)
    {
        id=i*(i-1)/2+j+3*m_N*(m_N+1)/2;
    }
    else if(j>i )//&& i<=m_N+1)
    {
        id=findid(0,m_N, m_N-i);
    }
    else if(i==m_N-1)
    {
        int n=m_N;
        int m=j+1;
        id=n*(n-1)/2+2*m_N*(m_N+1)/2+m;
    }
    else if(j==0)
        id=findid(1,m_N,i+2);        // id=0;
    
    return id;
    
}
