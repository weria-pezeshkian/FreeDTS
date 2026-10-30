#ifndef AFX_CreateMashBluePrint_H_INCLUDED_
#define AFX_CreateMashBluePrint_H_INCLUDED_

/*
 * @brief CreateMashBluePrint class
 *
 * Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 * This class creates a system using the provided inclusion definition in the input file and topology file.
 */

#include "inclusion.h"
#include "SimDef.h"

// Data structures for loading topology files into the RAM
struct Vertex_Map {    // Data structure for vertex map (not vertex object)
    double x, y, z;
    int id, domain;
    bool include;
};

struct Triangle_Map {    // Data structure for triangle map (not triangle object)
    int id, v1, v2, v3;
};

struct Inclusion_Map {    // Data structure for inclusion map 
    double x, y;
    int id, vid, tid;
};
struct VectorField_Map {    // Data structure for vector field map
    std::string data_line;
};

struct MeshBluePrint {    // Data structure for the mesh blueprint
    std::vector<Vertex_Map> bvertex;       // Vector of all vertices (only the blueprint, not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // Vector of all triangles (only the blueprint, not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // Vector of all inclusions (only the blueprint, not the object) in the mesh
    std::vector<int> excluded_id;
    Vec3D simbox;
    int number_vector_field;
    std::vector<VectorField_Map> bvectorfields;
};

class CreateMashBluePrint {
public:
    CreateMashBluePrint();
    ~CreateMashBluePrint();
    
    MeshBluePrint MashBluePrintFromInput_Top(const std::string& inputfilename, const std::string& topfilename); // Function to generate blueprint from provided topology file and input file
    
private:
    std::string m_InputFileName;  // Class variable for the name of the input file (dts file)
    std::string m_TopologyFileName; // Class variable for the name of the topology file

    bool ReadTopology(const std::string& file);  // Function to generate a mesh topology using the provided files. If the restart is on, then the topology will be generated from the restart file.
    void Read_Mult_QFile(const std::string& file);
    void Read_TSIFile(const std::string& topfile);
    bool GenerateIncFromInputfile(); // Function to generate some distribution of inclusions based on the input file. It does this only if the topology is from q files, since the tsi file format should have inclusions inside ...

    bool m_Healthy;   // To check if the input data are read correctly

    // Private members containing the mesh data (A stands for all)
    std::vector<Vertex_Map> m_VertexMap;       // Vector of all vertices (only the blueprint, not the object) in the mesh
    std::vector<Triangle_Map> m_TriangleMap;   // Vector of all triangles (only the blueprint, not the object) in the mesh
    std::vector<Inclusion_Map> m_InclusionMap; // Vector of all inclusions (only the blueprint, not the object) in the mesh
    std::vector<int> m_ExcludedID;
    Vec3D m_Box;
    std::vector<VectorField_Map> m_VectorFieldsMap;
    MeshBluePrint m_MeshBluePrint;
    int m_Number_of_VectorFields;
};

#endif
