
#include <time.h>
#include <iomanip>
#include "TSIFile.h"


TSIFile::TSIFile()
{
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);

}

TSIFile::~TSIFile()
{
    
}
void TSIFile::WriteTSI(std::string filename , MeshBluePrint blueprint)
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
MeshBluePrint TSIFile::ReadTSI(std::string filename)
{
    MeshBluePrint blueprint;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    Vec3D simbox;


    //======================
    int step,nver,ntr,ninc,id,tid,vid;
    double x,y,z;
    int v1,v2,v3,domain;
    std::string version,str;
    Nfunction f;
    std::ifstream input;
    input.open(filename.c_str());
    while (true)
    {
        input>>str;
        if(input.eof())
            break;
        
        if(str=="box")
        {
            getline(input,str);
            std::vector<std::string> S = f.split(str);
            if(S.size()<3)
            {
                std::cout<<"error ---> information of the box is not sufficent in the tsi file \n";
                break;
            }
            else
            {
                simbox(0) = f.String_to_Double(S.at(0));
                simbox(1) = f.String_to_Double(S.at(1));
                simbox(2) = f.String_to_Double(S.at(2));
            }
        }
        else if(str=="vertex")
        {
            input>>nver;
            double x,y,z;
            int id,domain;
            getline(input,str);

            for (int i=0;i<nver;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);

                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the vertex "<<i<<" is not sufficent in the tsi file \n";
                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                    Vertex_Map v;
                    id=f.String_to_Int(S.at(0));
                    x=f.String_to_Double(S.at(1));
                    y=f.String_to_Double(S.at(2));
                    z=f.String_to_Double(S.at(3));
                    v.id = id;
                    v.x = x;
                    v.y = y;
                    v.z = z;
                    v.domain = 0;
                    if(S.size()>4)
                    {
                        domain=f.String_to_Double(S.at(4));
                        v.domain = domain;
                    }
                    bvertex.push_back(v);
                }
            }

        }
        else if(str=="triangle")
        {
            input>>ntr;
            getline(input,str);

            int id,v1,v2,v3;

            for (int i=0;i<ntr;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the trinagles  "<<i<<" is not sufficent in the tsi file \n";
                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                    Triangle_Map t;
                    id=f.String_to_Int(S.at(0));
                    v1=f.String_to_Int(S.at(1));
                    v2=f.String_to_Int(S.at(2));
                    v3=f.String_to_Int(S.at(3));
                    t.id =id;
                    t.v1 = v1;
                    t.v2 = v2;
                    t.v3 = v3;
                    btriangle.push_back(t);
                    
                }
            }
        }
        else if(str=="inclusion")
        {
            input>>ninc;
            getline(input,str);

            int id,vid,type;
            double Dx, Dy;
            for (int i=0;i<ninc;i++)
            {
                getline(input,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<5)
                {
                    std::cout<<"error ---> information of the inclusion "<<i<<" is not sufficent in the tsi file \n";
                    std::cout<<str<<" \n";

                    break;
                }
                else
                {
                id=f.String_to_Int(S.at(0));
                type=f.String_to_Int(S.at(1));
                vid=f.String_to_Int(S.at(2));
                Dx=f.String_to_Double(S.at(3));
                Dy=f.String_to_Double(S.at(4));
                double norm = sqrt(Dx*Dx+Dy*Dy);
                    Dy=Dy/norm;
                    Dx=Dx/norm;
                    Inclusion_Map INC;
                    INC.x = Dx;
                    INC.y = Dy;
                    INC.id = id;
                    INC.tid = type;
                    INC.vid = vid;
                    binclusion.push_back(INC);
                }

            }

        }
        else
        {
            getline(input,str);

        }
    }
    blueprint.bvertex = bvertex;
    blueprint.btriangle = btriangle;
    blueprint.binclusion = binclusion;
    blueprint.simbox = simbox;
    input.close();
    return blueprint;
}


