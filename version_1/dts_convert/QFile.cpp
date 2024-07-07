
#include <time.h>
#include <iomanip>
#include "QFile.h"


QFile::QFile()
{
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);

}

QFile::~QFile()
{
    
}
void QFile::Write(std::string filename , MeshBluePrint blueprint)
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
MeshBluePrint QFile::Read(std::string filename)
{
    MeshBluePrint blueprint;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    Vec3D Box;


    std::string str;
    Nfunction f;
    int vid = 0;
    int tid = 0;
    int id;

    std::ifstream Qs;
    Qs.open(filename.c_str());
    
    // first line is the box size and it should only contain 3 numbers;
    getline(Qs,str);
    std::vector<std::string> b = f.split(str);
    if(b.size()>3)
    {
        std::cout<<"---> Error: box information in the file "<<filename<<" is not correct "<<std::endl;
        exit(0);
    }
    // The final box size will be the largest box in all the q files
        Box(0)=f.String_to_Double(b[0]);
        Box(1)=f.String_to_Double(b[1]);
        Box(2)=f.String_to_Double(b[2]);

    // reading the number of the vertices in this file
    getline(Qs,str);
    b.clear();
    b = f.split(str);
    if(b.size()>1)
    {
        std::cout<<"----> Error: number of vertices in the file "<<filename<<" is not correct "<<std::endl;
        exit(0);
    }
    int NV = f.String_to_Int(b[0]);
    for (int i=0;i<NV;i++)
    {
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>5 || b.size()<4)
        {
            std::cout<<"----> Error: Line "<<i+2<<", info of a vertex in the file "<<filename<<" is not correct.  "<<std::endl;
            exit(0);
        }
        Vertex_Map v;
        v.id=vid;
        v.x=f.String_to_Double(b[1]);
        v.y=f.String_to_Double(b[2]);
        v.z=f.String_to_Double(b[3]);
        if(b.size()==5)
            v.domain=f.String_to_Int(b[4]);
        else
            v.domain=0;
        bvertex.push_back(v);
        vid++;
    }

    getline(Qs,str);
    b.clear();
    b = f.split(str);
    if(b.size()>1)
    {
        std::cout<<"----> Error: number of triangle in the file "<<filename<<" is not correct "<<str<<std::endl;
        exit(0);
    }
    int nt=f.String_to_Int(b[0]);
    for (int i=0;i<nt;i++)
    {
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>5 || b.size()<4)
        {
            std::cout<<"----> Error: Line "<<i+2<<", info of a triangle in the file "<<filename<<" is not correct.  "<<std::endl;
            exit(0);
        }
        Triangle_Map t;
        t.id=tid;
        t.v1=f.String_to_Int(b[1]);
        t.v2=f.String_to_Int(b[2]);
        t.v3=f.String_to_Int(b[3]);
        btriangle.push_back(t);
        tid++;
    }
    Qs.close();
    
    blueprint.bvertex = bvertex;
    blueprint.btriangle = btriangle;
    blueprint.binclusion = binclusion;
    blueprint.simbox = Box;
    
    
    return blueprint;
}


