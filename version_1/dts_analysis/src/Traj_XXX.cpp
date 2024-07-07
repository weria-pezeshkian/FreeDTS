
#include <time.h>
#include <iomanip>
#include "Traj_XXX.h"


Traj_XXX::Traj_XXX(Vec3D *pBox, std::string tsiFolder_name, std::string wORr)
{
 m_pBox=pBox;
m_Condition=true;
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);
    
    m_tsiFolder_name = tsiFolder_name;
  if(wORr=="w")
  {
    const int dir_err = system(("mkdir -p "+tsiFolder_name).c_str());
    if (-1 == dir_err)
    {
        std::cout<<"---> Error: while creating directory  "<<tsiFolder_name<<"\n";
        exit(1);
    }
  }
}

Traj_XXX::~Traj_XXX()
{
    
}
void Traj_XXX::WriteTSI(int step ,  std::string filename , std::vector< vertex* > pver, std::vector< triangle* > ptriangle,  std::vector< inclusion* > pinc)
{
    FILE * output;
    output = fopen((m_tsiFolder_name+"/"+filename).c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
//------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(*m_pBox)(0),(*m_pBox)(1),(*m_pBox)(2));

    const char* ver="vertex";
    int size=pver.size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<vertex *>::iterator it = pver.begin() ; it != pver.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetVID(),(*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
    
    const char* tri="triangle";
    size = ptriangle.size();
    fprintf(output,"%s%20d\n",tri,size);
    for (std::vector<triangle *>::iterator it = ptriangle.begin() ; it != ptriangle.end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",(*it)->GetTriID(),((*it)->GetV1())->GetVID(),((*it)->GetV2())->GetVID(),((*it)->GetV3())->GetVID());
    
    
    const char* inc="inclusion";
    size = pinc.size();
    fprintf(output,"%s%20d\n",inc,size);
    format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<inclusion *>::iterator it = pinc.begin() ; it != pinc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),(*it)->GetInclusionTypeID(),((*it)->Getvertex())->GetVID(),((*it)->GetLDirection())(0),((*it)->GetLDirection())(1));
 
    
    fclose(output);
    
}
MeshBluePrint Traj_XXX::ReadTSI(std::string filename,std::vector<InclusionType> bINCtype, bool *isok)
{
    Nfunction f;
    MeshBluePrint blueprint;
    
if (f.FileExist(filename))
{
    std::ifstream input;
    input.open(filename.c_str());
    std::string version = "0.0";
    std::string str;


    while (true)
    {
        getline(input,str);
        if(input.eof())
            break;
        std::vector<std::string> S = f.split(str);
    
        if(S.size()>1)
        if(S.at(0)=="version")
        {
            version = S.at(1);
        }
    }
    input.close();
    if(version=="0.0")
    {
        blueprint = ReadTSI1(filename);
    }
    else
    {
        blueprint = ReadTSI2(filename, bINCtype);

    }
    *isok = true;
}
else
{
    *isok = false;
}
    return blueprint;
}
    
MeshBluePrint Traj_XXX::ReadTSI2(std::string filename, std::vector<InclusionType> bINCtype)
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
                    m_Condition =false;
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
                    m_Condition =false;

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
                    m_Condition =false;

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
    blueprint.binctype = bINCtype;
    blueprint.simbox = simbox;
    input.close();
    return blueprint;
}
MeshBluePrint Traj_XXX::ReadTSI1(std::string filename)
{
    MeshBluePrint blueprint;
    //=======================
    std::cout<<"----> Error: your tsi file is from an old version, we do not support that version anymore "<<std::endl;
    return blueprint;
}


