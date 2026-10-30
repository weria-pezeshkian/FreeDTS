#ifndef TriangularPrismBuilder_H
#define TriangularPrismBuilder_H
#include <array>
#include "SimDef.h"
#include "triangle.h"
#include "Vec3D.h"
#include "Tensor2.h"


/*
 * -----------------------------------------------------------------------------
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Date: June 2026
 *
 * Class: TriangularPrismBuilder
 *
 * Description:
 * Constructs and analyzes all possible triangular-prism-like connections
 * between two triangular faces.
 *
 * Given two triangles, the builder generates the possible vertex
 * correspondences and side-face triangulations that connect them into a
 * closed surface. The generated configurations can be used for geometric,
 * topological, or mesh-construction purposes.
 *
 * The class stores references to the two input triangles and the simulation
 * box required for geometric operations involving periodic boundary
 * conditions.
 *
 * Notes:
 * - Multiple valid prism configurations may exist for the same triangle pair.
 * - The class does not own the input triangles or simulation box.
 *
 * -----------------------------------------------------------------------------
 * License:
 * Copyright (c) Weria Pezeshkian, 2026. All rights reserved.
 * -----------------------------------------------------------------------------
 */
using PairMap = std::array<int, 2>; 
inline PairMap make_PairMap(int a, int b) {
    return {a,b};
} 
// At namespace scope
struct PairMapHash {
    size_t operator()(const PairMap& p) const {
        return std::hash<int>()(p[0]) ^ (std::hash<int>()(p[1]) << 1);
    }
};
using triple = std::array<int, 3>;
inline triple make_triple(int a, int b, int c) {
    return {a,b,c};
}
struct TriangularPrism {
    std::vector<triple> VTriples;
    std::vector<Vec3D> VNormals;
};
class TriangularPrismBuilder {
public:

    /*
     * Construct a prism builder from two triangles.
     *
     * Parameters:
     *   id   - Unique identifier.
     *   t1   - First triangle.
     *   t2   - Second triangle.
     *   pbox - Simulation box used for geometric calculations.
     *   max distance 
     */
    TriangularPrismBuilder(Vec3D* pbox, double &maxl2, double &minagle, std::string &topfilename);

    ~TriangularPrismBuilder();


public:

        std::vector<TriangularPrism> GeneratePossibleTopology(triangle* t1,  triangle* t2);

private:
void GenerateDistanceMap(triangle* t1,  triangle* t2);  
bool CheckMap(TriangularPrism & tp);
bool CheckDistanceEdge(int i, int j);    
bool MakePrismMaps();  
bool ReorderTriple(triple& t);
Vec3D TripleNormal(triple& t);
bool TriplesAreNeighbour(triple& t1, triple& t2);
bool UpdatePrismNormals(TriangularPrism & tp);
bool CheckPrismNormals(TriangularPrism & tp);
bool MakeDefaultPrismMaps();
TriangularPrism ReadPrismMap(std::istream& input);
bool CheckPrismMapOrientation(TriangularPrism & tp); // Checks whether the triangular prism face map can be made with consistent orientation.
private:

    //--------------------------------------------------------------------------
    // Data members
    //--------------------------------------------------------------------------
     Vec3D* m_pBox; // Ownership is not required.
    double m_Dist[3][3];
    std::vector<TriangularPrism> m_TopologyMaps;
    double &m_MaxLength2;
    // just to keep as global
    std::vector<vertex *> m_Vertices;
    double &m_MinAngle;
    
     std::string &m_TopologyFilename;
    

};

#endif
