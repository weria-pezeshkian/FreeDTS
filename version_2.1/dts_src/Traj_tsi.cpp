
#include "State.h"
#include "Vec3D.h"
#include "triangle.h"
#include "vertex.h"
#include "Traj_tsi.h"
#include "Nfunction.h"
Traj_tsi::Traj_tsi(State *pstate){

    m_Period = 1000;
    m_Folder_name = "TrajTSI";
    m_pState = pstate;
    m_Precision = TSI_Precisions;
}
Traj_tsi::Traj_tsi(State *pstate, int period, std::string tsiFolder_name){
    
    m_pState = pstate;
    m_Period = period;
    m_Folder_name = tsiFolder_name;
    m_Precision = TSI_Precisions;
    
}

Traj_tsi::~Traj_tsi(){
    
}
bool Traj_tsi::OpenFolder() {
    // Use Nfunction class for centralized file handling functions
    // waiting for use of <filesystem> in c++17; but lets not higher up the required version of c++ for other users
    return  Nfunction::OpenFolder(m_Folder_name);
}
void Traj_tsi::WriteAFrame(int step){
 
    // If the period is not the right step, ignore the call
    if (m_Period == 0 || step % m_Period != 0)
        return;
    
    int id = step/m_Period;
    
    std::string filename = m_pState->GetRunTag() + Nfunction::D2S(id)+ "." + TSIExt;
    
    WriteAFrame(filename);
    
    return;
}
void Traj_tsi::WriteAFrame(std::string filename){

    // Check the extension of the filename and add ".tsi" if needed
    if (Nfunction::SubstringFromRight(filename,'.') != TSIExt)
        filename =filename + "." + TSIExt;
    
    // Retrieve active vertices, triangles, and inclusions from the mesh
     std::vector<vertex*> pver = m_pState->GetMesh()->GetActiveV();
     std::vector<triangle*> ptriangle = m_pState->GetMesh()->GetActiveT();
     std::vector<inclusion*> pinc = m_pState->GetMesh()->GetInclusion();
    Vec3D *pBox = m_pState->GetMesh()->GetBox();
    // Open the file for writing
    FILE* output = fopen((m_Folder_name + "/" + filename).c_str(), "w");
    if (!output) {
        std::cerr << "Failed to open file for writing: " << m_Folder_name << "/" << filename << std::endl;
        return;
    }

    std::string format = "%"+m_Precision+"lf%"+m_Precision+"lf%"+m_Precision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
//------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(*pBox)(0),(*pBox)(1),(*pBox)(2));

    const char* ver="vertex";
    int size=pver.size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_Precision+"lf%"+m_Precision+"lf%"+m_Precision+"lf\n";
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
    format = "%10d%10d%10d%"+m_Precision+"lf%"+m_Precision+"lf\n";
    for (std::vector<inclusion *>::iterator it = pinc.begin() ; it != pinc.end(); ++it)
        fprintf(output,format.c_str(),(*it)->GetID(),(*it)->GetInclusionType()->ITid,((*it)->Getvertex())->GetVID(),((*it)->GetLDirection())(0),((*it)->GetLDirection())(1));
 
    
    size = m_pState->GetMesh()->GetNoVFPerVertex();
    if(size != 0){
        const char* vfname = "vector_fields";
        fprintf(output, "%s%20d\n", vfname, size);
        for (std::vector<vertex*>::iterator it = pver.begin(); it != pver.end(); ++it) {
            std::string s_vf = (*it)->GetVectorFieldsStream();
            fprintf(output, "%s\n", s_vf.c_str());
        }
    }
    
    fclose(output);
    
}
/*
MeshBluePrint Traj_tsi::ReadAFrame(std::string filename, bool &isok)
{
    Nfunction f;
    MeshBluePrint blueprint;
    
    
    std::vector<InclusionType> bINCtype;
    std::cout<<" this function is now complete \n";
    
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
        std::cout<<"----> Error: your tsi file is from an old version, we do not support that version anymore "<<std::endl;
    }
    else
    {
        blueprint = ReadTSI2(filename, bINCtype);

    }
    isok = true;
}
else
{
    isok = false;
}
    return blueprint;
}
    
MeshBluePrint Traj_tsi::ReadTSI2(std::string filename, std::vector<InclusionType> bINCtype)
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

*/
std::string Traj_tsi::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+m_Folder_name+" "+Nfunction::D2S(m_Period);
    return state;
}
