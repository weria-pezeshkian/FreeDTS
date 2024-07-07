

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "CreateMashBluePrint.h"
#include "Nfunction.h"
#include "RNG.h"
CreateMashBluePrint::CreateMashBluePrint()
{
}
MeshBluePrint CreateMashBluePrint::MashBluePrintFromInput_Top(std::string inputfilename, std::string topfilename)
{
// initializing the box size
    m_Box(0) = 1;
    m_Box(1) = 1;
    m_Box(2) = 1;

#if TEST_MODE == Enabled
    std::cout<<"----> We have reached the CreateMashBluePrint Class -- "<<std::endl;
#endif
    m_InputFileName = inputfilename;
    m_TopologyFileName = topfilename;
    
    //==== We read the inclusions type from the input file, we create inclusion types
    ReadInclusionType(m_InputFileName);
    
#if TEST_MODE == Enabled
    std::cout<<"----> Inclusion types has been read from the input file -- "<<std::endl;
#endif
    //== We read the topology file, this could be in a several q file, or a tsi file.
#if TEST_MODE == Enabled
    std::cout<<"----> atempted to read topology of the system from the createmashblueprintclass -- "<<std::endl;
#endif
    ReadTopology(topfilename);
#if TEST_MODE == Enabled
    std::cout<<"----> topology of the system was read from tsi or q in the createmashblueprintclass -- "<<std::endl;
#endif
    m_MeshBluePrint.bvertex = m_VertexMap;
    m_MeshBluePrint.btriangle = m_TriangleMap;
    m_MeshBluePrint.binclusion = m_InclusionMap;
    m_MeshBluePrint.binctype = m_AllInclusionType;
    m_MeshBluePrint.simbox = m_Box;
    
    return m_MeshBluePrint;
}
CreateMashBluePrint::~CreateMashBluePrint()
{
    
}
void CreateMashBluePrint::ReadTopology(std::string file)
{

        std::string ext = file.substr(file.find_last_of(".") + 1);
        if(ext=="tsi")
        {
            Read_TSIFile(file);

        }
        else if(ext=="top")
        {
            Read_Mult_QFile(file);
            GenerateIncFromInputfile();
        }


}
void CreateMashBluePrint::Read_TSIFile(std::string tsifile)
{
    Nfunction f;
    std::ifstream tsi;
    tsi.open(tsifile.c_str());
    std::string str;
    int nver,ntr,ninc;
    while (true)
    {
        tsi>>str;
        if(tsi.eof())
            break;
        if(str=="version")
        {
            getline(tsi,str);
        }
        else if(str=="box")
        {
            getline(tsi,str);
            std::vector<std::string> S = f.split(str);
            if(S.size()<3)
            {
                std::cout<<"---> Error, information of the box is not sufficent in the tsi file \n";
                exit(0);
            }
            else
            {
                m_Box(0) = f.String_to_Double(S.at(0));
                m_Box(1) = f.String_to_Double(S.at(1));
                m_Box(2) = f.String_to_Double(S.at(2));
            }
#if DEBUG_MODE == Enabled
            std::cout<<"----> box was read  "<<std::endl;
            std::cout<<m_Box(0)<<"  "<<m_Box(1)<<"  "<<m_Box(2)<<"  "<<str<<std::endl;
#endif
        }
        else if(str=="vertex")
        {
            tsi>>nver;
            getline(tsi,str);
            
            for (int i=0;i<nver;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                
                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the vertex "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Vertex_Map v;
                    v.id=i;
                    v.x=f.String_to_Double(S[1]);
                    v.y=f.String_to_Double(S[2]);
                    v.z=f.String_to_Double(S[3]);
                    if(S.size()>4)
                        v.domain=f.String_to_Int(S.at(4));
                    else
                        v.domain=0;
                    m_VertexMap.push_back(v);
                }
            }
            
        }
        else if(str=="triangle")
        {
            tsi>>ntr;
            getline(tsi,str);
            for (int i=0;i<ntr;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<4)
                {
                    std::cout<<"error ---> information of the trinagles  "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Triangle_Map t;
                    t.id=f.String_to_Int(S.at(0));
                    t.v1=f.String_to_Int(S.at(1));
                    t.v2=f.String_to_Int(S.at(2));
                    t.v3=f.String_to_Int(S.at(3));
                    m_TriangleMap.push_back(t);
                    
                }
            }
        }
        else if(str=="inclusion")
        {
            tsi>>ninc;
            getline(tsi,str);
            for (int i=0;i<ninc;i++)
            {
                getline(tsi,str);
                std::vector<std::string> S = f.split(str);
                if(S.size()<5)
                {
                    std::cout<<"error ---> information of the inclusion "<<i<<" is not sufficent in the tsi file \n";
                    exit(0);
                }
                else
                {
                    Inclusion_Map inc;
                    inc.id=f.String_to_Int(S.at(0));
                    inc.tid=f.String_to_Int(S.at(1));
                    inc.vid=f.String_to_Int(S.at(2));
                    double x=f.String_to_Double(S.at(3));
                    double y=f.String_to_Double(S.at(4));
                    double norm = sqrt(x*x+y*y);
                    x=x/norm;
                    y=y/norm;
                    inc.x = x; inc.y = y;
                    m_InclusionMap.push_back(inc);
                }
            }
        }
        else if(str=="exclusion")
        {
            break;
        }
        else
        {
            std::cout<<"error ---> "<<str<<" is unidentified key word for tsi file \n";
            exit(0);
        }
    }
}
void CreateMashBluePrint::Read_Mult_QFile(std::string topfile)
{
    Nfunction f;
    //== read the top file and store all the q files with the group name.
    std::vector<std::string> qfiles;
    std::vector<int> groupid;
    std::string str;
    int id;
    std::ifstream top;
    top.open(topfile.c_str());
    while (true)
    {
        top>>str>>id;
        if(top.eof())
            break;
        qfiles.push_back(str);
        groupid.push_back(id);
#if TEST_MODE == Enabled
        std::cout<<"----> q file in the top file: "<<str<<std::endl;
#endif
        getline(top,str);
    }
    top.close();
    // read each q file
    double Lx,Ly,Lz;
    int vid = 0;
    int tid = 0;
    for (int fi=0;fi<qfiles.size();fi++)
    {
        std::ifstream Qs;
        Qs.open((qfiles.at(fi)).c_str());
        
        // first line is the box size and it should only contain 3 numbers;
        getline(Qs,str);
        std::vector<std::string> b = f.split(str);
        if(b.size()>3)
        {
            std::cout<<"---> Error: box information in the file "<<qfiles.at(fi)<<" is not correct "<<std::endl;
            exit(0);
        }
        // The final box size will be the largest box in all the q files
        if(m_Box(0)<f.String_to_Double(b[0]))
            m_Box(0)=f.String_to_Double(b[0]);
        if(m_Box(1)<f.String_to_Double(b[1]))
            m_Box(1)=f.String_to_Double(b[1]);
        if(m_Box(2)<f.String_to_Double(b[2]))
            m_Box(2)=f.String_to_Double(b[2]);
#if DEBUG_MODE == Enabled
        std::cout<<"----> box was read  "<<std::endl;
        std::cout<<m_Box(0)<<"  "<<m_Box(1)<<"  "<<m_Box(2)<<"  "<<str<<std::endl;
#endif
        // reading the number of the vertices in this file
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>1)
        {
            std::cout<<"----> Error: number of vertices in the file "<<qfiles.at(fi)<<" is not correct "<<std::endl;
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
                std::cout<<"----> Error: Line "<<i+2<<", info of a vertex in the file "<<qfiles.at(fi)<<" is not correct.  "<<std::endl;
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
            m_VertexMap.push_back(v);
            vid++;
        }
#if DEBUG_MODE == Enabled
        std::cout<<"----> vertex section was read  "<<std::endl;
#endif
        getline(Qs,str);
        b.clear();
        b = f.split(str);
        if(b.size()>1)
        {
            std::cout<<"----> Error: number of triangle in the file "<<qfiles.at(fi)<<" is not correct "<<str<<std::endl;
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
                std::cout<<"----> Error: Line "<<i+2<<", info of a triangle in the file "<<qfiles.at(fi)<<" is not correct.  "<<std::endl;
                exit(0);
            }
            Triangle_Map t;
            t.id=tid;
            t.v1=f.String_to_Int(b[1]);
            t.v2=f.String_to_Int(b[2]);
            t.v3=f.String_to_Int(b[3]);
            m_TriangleMap.push_back(t);
            tid++;
        }
        Qs.close();
    }
    std::cout<<"trinagle is read "<<"\n";

}
void CreateMashBluePrint::ReadInclusionType(std::string file)
{
    Nfunction f;
    // enforcing correct file extension:  check simdef file for value of InExt
    std::string ext = file.substr(file.find_last_of(".") + 1);
    if(ext!=InExt)
    file = file + "." + InExt;
    if (f.FileExist(file)!=true)
    {
        std::cout<<"----> Error: the input file with the name "<<file<< " does not exist "<<std::endl;
        m_Healthy =false;
        exit(0);
    }
    std::ifstream input;
    input.open(file.c_str());
    std::string firstword,rest,str1,str2,TypeNames;
    int N,TypeID, NoType;
    double Kappa,KappaG,KappaP,KappaL,C0,C0P,C0N;
    InclusionType EmptyIncType;
    {
        EmptyIncType.ITName = "noinc";   // type name
        EmptyIncType.ITid = 0;       // Type ID
        EmptyIncType.ITN = 0;      // inplane symmetry
        EmptyIncType.ITk =0;     // Type Kappa
        EmptyIncType.ITkg =0;  // kappaG
        EmptyIncType.ITk1 =0;  // K_||
        EmptyIncType.ITk2 =0;  // K_norm
        EmptyIncType.ITc0 =0;  // curvature
        EmptyIncType.ITc1 =0;    // direction curvature 1
        EmptyIncType.ITc2 =0;    // direction curvature 2
    }
    m_AllInclusionType.push_back(EmptyIncType); // we will make a default inclusion type
    while (true)
    {
        input>>firstword;
        if(input.eof() || firstword == "INCLUSION")
            break;
    }

            input>>str1>>NoType>>str2;
            getline(input,rest);
            getline(input,rest); // takes a text line of SRotation Type   K  KG KP KL C0 C0P C0L no restriction for it yet

#if TEST_MODE == Enabled
    std::cout<<" --- > This line should be as below "<<rest<<std::endl;
    std::cout<<" --- >                              SRotation Type   K  KG KP KL C0 C0P C0L "<<std::endl;
#endif
    
    if(str1=="Define" || str1=="define"){
        for(int i=0;i<NoType;i++)
        {
        input>>N>>TypeNames>>Kappa>>KappaG>>KappaP>>KappaL>>C0>>C0P>>C0N;
        InclusionType IncType;
            {
                IncType.ITName = TypeNames;   // type name
                IncType.ITid = i+1;       // Type ID
                IncType.ITN = N;      // inplane symmetry
                IncType.ITk = Kappa;     // Type Kappa
                IncType.ITkg = KappaG;  // kappaG
                IncType.ITk1 = KappaP;  // K_||
                IncType.ITk2 = KappaL;  // K_norm
                IncType.ITc0 = C0;  // curvature
                IncType.ITc1 = C0P;    // direction curvature 1
                IncType.ITc2 = C0N;    // direction curvature 2
            }
        m_AllInclusionType.push_back(IncType);
        }
    }
    input.close();
        WriteCreateMashBluePrintLog();

}
void CreateMashBluePrint::GenerateIncFromInputfile()
{
    std::string file = m_InputFileName;
    Nfunction f;
    std::string str;
    // enforcing correct file extension:  check simdef file for value of InExt
    std::string ext = file.substr(file.find_last_of(".") + 1);
    if(ext!=InExt)
        file = file + "." + InExt;
    if (f.FileExist(file)!=true)
    {
        std::cout<<"----> Error: the input file with the name "<<file<< " does not exist "<<std::endl;
        m_Healthy =false;
        exit(0);
    }
    std::ifstream input;
    input.open(file.c_str());
    std::vector<int>  Type;
    std::vector<double> Density;
     while (true)
    {
        input>>str;
        if(input.eof())
            break;
        if(str=="GenerateInclusions")
        {
            input>>str>>str;
            if(str=="Random")
            {
                getline(input,str);
                std::string type, den;
                getline(input,type);
                getline(input,den);
                std::vector <std::string> T = f.split(type);
                std::vector <std::string> D = f.split(den);
                double totdensity = 0;
                for (int i=1;i<T.size();i++)
                {
                    Type.push_back(f.String_to_Int(T.at(i)));
                    Density.push_back(f.String_to_Double(D.at(i)));
                    totdensity+=f.String_to_Double(D.at(i));
                }
                if(T.size()!=D.size() || T.at(0)!="TypeID" || D.at(0)!="Density" || totdensity>1)
                {
                    std::cout<<"----> Error: info in inclusion type and density are not correct "<<std::endl;
#if TEST_MODE == Enabled
                    std::cout<<T.size()<<"  "<<D.size()<<"  "<<T.at(0)<<"   "<<D.at(0)<<std::endl;
#endif
                    exit(0);
                }
                break;
            }
            else
            {
                std::cout<<"----> Error: "<<str<<" type of the distribution is not defined yet "<<std::endl;
                exit(0);
            }
        
        }
        else
        {
            getline(input,str);
        }
    }//while(true)
    int *V = new int[m_VertexMap.size()];  // A variable
    for (int i=0;i<m_VertexMap.size();i++)
        V[i]=0;
    double incid = 0;

    RNG rng(3434);
                for (int j=0;j<Density.size();j++)
                {
                    double Createdinc = 0;
                    int Ninc = int(double(m_VertexMap.size())*Density.at(j));
                    int iter = 0;
                    while(true)
                    {
                        iter++;
                        if(Createdinc==Ninc || iter>100*(m_VertexMap.size()))
                            break;
                        int n = rng.UniformRNG(m_VertexMap.size());
                        if(V[n]==0)
                        {
                            double x=1-2*(rng.UniformRNG(1));
                            double y=1-2*(rng.UniformRNG(1));
                            double r=sqrt(x*x+y*y);
                            x=x/r;
                            y=y/r;
                            V[n]=1;
                            Inclusion_Map inc;
                            inc.id=incid;
                            inc.tid=Type[j];
                            inc.vid=n;
                            inc.x = x; inc.y = y;
                            m_InclusionMap.push_back(inc);
                            incid++;
                            Createdinc++;
                        }
                    }
                }
    delete [] V;
}
void CreateMashBluePrint::WriteCreateMashBluePrintLog()
{
    std::ofstream log;
    log.open("statelog.log",std::fstream::app);
    log<<"INCLUSION"<<std::endl;
    log<<"Define "<<m_AllInclusionType.size()<<" Inclusions "<<std::endl;
    log<<" SRotation Type   K  KG KP KL C0 C0P C0L "<<std::endl;
    for (std::vector<InclusionType>::iterator it = m_AllInclusionType.begin() ; it != m_AllInclusionType.end(); ++it)
        log<<it->ITN<<"  "<<it->ITName<<"  "<<it->ITk<<"  "<<"  "<<it->ITkg<<"  "<<"  "<<it->ITk1<<"  "<<"  "<<it->ITk2<<"  "<<"  "<<it->ITc0<<"  "<<"  "<<it->ITc1<<"  "<<"  "<<it->ITc2<<"  "<<std::endl;
    log.close();
}
