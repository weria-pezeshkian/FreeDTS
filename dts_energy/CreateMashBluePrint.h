#if !defined(AFX_CreateMashBluePrint_H_INCLUDED_)
#define AFX_CreateMashBluePrint_H_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class creats a system using the provided inclusion definition in the input file and topology file or restart file or etc.
 */
#include "inclusion.h"
#include "SimDef.h"
// this three data structure is for loading topology files into the RAM 
struct Vertex_Map {    // data structure for vertex map (not vertex object)
    double x,y,z;
    int id,domain;
};
struct Triangle_Map {    // data structure for triangle map (not triangle object)
    int id,v1,v2,v3;
};
struct Inclusion_Map {    // data structure for inclusion map (not triangle object)
    double x,y;
    int id,vid,tid;
};
struct MeshBluePrint {    // data structure for the mesh blue print
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector <InclusionType> binctype;  // a vector containing all inclsuion type and a default one
    Vec3D simbox;
};
class CreateMashBluePrint
{
public:
	CreateMashBluePrint();
	 ~CreateMashBluePrint();
    
public:
    
    std::vector <InclusionType> GetAllInclusionType()    const {return m_AllInclusionType;}
public:
    MeshBluePrint MashBluePrintFromInput_Top(std::string inputfilename, std::string topfilename); // A function to generate bluee print from provided topology file and input file //
    void ReadInclusionType(std::string file);  // a function to read the inclsuion section of the input file.
private:
    std::string m_InputFileName;  // a class  variable for the name of the input file (dts file)
    std::string m_TopologyFileName; // a class  variable for the name of the topology file
private:
    void WriteCreateMashBluePrintLog();  // writing a log file about what has been read and what is out
    void ReadTopology(std::string file);  // a function to generate a mesh topology using the provided files. If the restart is on, then the topology will be generated from the restart file.
    void Read_Mult_QFile(std::string);
    void Read_TSIFile(std::string topfile);
    void GenerateIncFromInputfile(); // this function generate some distribution of inclsuions based on the input file. It do this only if the topology is from q files, since the tsi file format should have inclusion inside ...

private:
    bool m_Healthy;   // To check if the input data are read correctly

  // private memebrs containing the mesh data (A) stands for all
private:
    std::vector<Vertex_Map> m_VertexMap;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> m_TriangleMap;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> m_InclusionMap; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector <InclusionType> m_AllInclusionType;  // a vector containing all inclsuion type and a default one
    Vec3D m_Box;
private:
    MeshBluePrint m_MeshBluePrint;

};

#endif
