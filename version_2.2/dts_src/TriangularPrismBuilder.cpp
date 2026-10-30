#include <fstream>
#include <sstream>
#include "TriangularPrismBuilder.h"
#include "vertex.h"
#include "MESH.h"

// Constructor with initialization of TriangularPrismBuilder parameters
TriangularPrismBuilder::TriangularPrismBuilder(Vec3D *pbox, double &maxl2, double &minagle, std::string &topfilename)
    : m_pBox(pbox), m_MaxLength2(maxl2), m_MinAngle(minagle), m_TopologyFilename(topfilename) {

        if(!MakePrismMaps()){
            std::cout<<" error--> 220202 \n";
        }
}
// Destructor
TriangularPrismBuilder::~TriangularPrismBuilder() = default;
std::vector<TriangularPrism> TriangularPrismBuilder::GeneratePossibleTopology(triangle* t1,  triangle* t2){
    
std::vector<TriangularPrism> ValidTopology;

            GenerateDistanceMap(t1 , t2); // as this is constant for the two trinagle, we only calculate it once
 
            for (TriangularPrism& prism : m_TopologyMaps) {
                    if (CheckMap(prism)){
                            ValidTopology.push_back(prism);
                    }
            }
    
    return ValidTopology;
}

bool TriangularPrismBuilder::CheckMap(TriangularPrism &tprism) {
    
    // check distances
    for (const auto& tri : tprism.VTriples)
    {
        for (int k = 0; k < 3; k++)
        {
            int i = tri[k];
            int j = tri[(k + 1) % 3];

            if (!CheckDistanceEdge(i, j))
                return false;
        }
    }
    
    // check faces
      if(!UpdatePrismNormals(tprism)){
             return false;
      }
    
     if(!CheckPrismNormals(tprism)){
             return false;
      }

    return true;
}
bool TriangularPrismBuilder::CheckPrismNormals(TriangularPrism &tp)
{
     auto& triples = tp.VTriples;
     auto& normals = tp.VNormals;
    const size_t n = triples.size();

    // -----------------------------
    // 1. Check neighbouring triples
    // -----------------------------
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < n; j++)
        {
            if (!TriplesAreNeighbour(triples[i], triples[j]))
                continue;

            const double angle = Vec3D::dot(normals[i], normals[j]);

            if (angle < m_MinAngle){
                return false;
            }
        }
    }

    // -----------------------------------------
    // 2. Check against original mesh triangles
    // -----------------------------------------
    for (size_t i = 0; i < n; i++)
    {
        const triple& tri = triples[i];

        vertex* v1 = m_Vertices[tri[0]];
        vertex* v2 = m_Vertices[tri[1]];

        links* mylink = nullptr;

        for (auto* l : v1->GetVLinkList())
        {
            if (l->GetV2() == v2)
            {
                mylink = l;
                break;
            }
        }

        if (!mylink || !(mylink->GetMirrorFlag())){
            std::cout <<"---> error (id 4647393), this should not happen \n";
            return false; // safety guard
        }

        const double angle  = Vec3D::dot(normals[i], mylink->GetMirrorLink()->GetTriangle()->GetNormalVector());

        if (angle < m_MinAngle){
            return false;
        }
    }

    return true;
}
bool TriangularPrismBuilder::TriplesAreNeighbour(triple& t1, triple& t2)
{
    int matchCount = 0;

    // check each element of t1 against t2
    if (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2]) matchCount++;
    if (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2]) matchCount++;
    if (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2]) matchCount++;

    return matchCount == 2;
}
bool TriangularPrismBuilder::CheckDistanceEdge(int i, int j)
{
    // same-triangle edges (optional skip)
    if ((i < 3 && j < 3) || (i >= 3 && j >= 3))
        return true;

    int ia = (i < 3) ? i : i - 3;
    int jb = (j < 3) ? j : j - 3;

    return (m_Dist[ia][jb] <= m_MaxLength2);
}
void TriangularPrismBuilder::GenerateDistanceMap(triangle* t1,  triangle* t2){
    
        // ---- vertex sharing rejection ----
     vertex * v1 = t1->GetV1();
     vertex * v2 = t1->GetV2();
     vertex * v3 = t1->GetV3();

     vertex * u1 = t2->GetV1();
     vertex * u2 = t2->GetV2();
     vertex * u3 = t2->GetV3();
     
    m_Vertices.clear();
    m_Vertices.insert(m_Vertices.end(),{v1, v2, v3, u1, u2, u3});

    double dist_11 = MESH::SquareDistanceBetweenTwoVertices(v1, u1, *m_pBox);
    double dist_12 = MESH::SquareDistanceBetweenTwoVertices(v1, u2, *m_pBox);
    double dist_13 = MESH::SquareDistanceBetweenTwoVertices(v1, u3, *m_pBox);
    double dist_21 = MESH::SquareDistanceBetweenTwoVertices(v2, u1, *m_pBox);
    double dist_22 = MESH::SquareDistanceBetweenTwoVertices(v2, u2, *m_pBox);
    double dist_23 = MESH::SquareDistanceBetweenTwoVertices(v2, u3, *m_pBox);
    double dist_31 = MESH::SquareDistanceBetweenTwoVertices(v3, u1, *m_pBox);
    double dist_32 = MESH::SquareDistanceBetweenTwoVertices(v3, u2, *m_pBox);
    double dist_33 = MESH::SquareDistanceBetweenTwoVertices(v3, u3, *m_pBox);

    m_Dist[0][0] = dist_11;
    m_Dist[0][1] = dist_12;
    m_Dist[0][2] = dist_13;
    m_Dist[1][0] = dist_21;
    m_Dist[1][1] = dist_22;
    m_Dist[1][2] = dist_23;
    m_Dist[2][0] = dist_31;
    m_Dist[2][1] = dist_32;
    m_Dist[2][2] = dist_33;
    
    
    return;
}
bool TriangularPrismBuilder::ReorderTriple(triple& t) {
// Reorders a triple of values in the range [0..5] according to two cyclic orientation groups:
//
//   Group A: {0,1,2} with cycle 0 → 1 → 2 → 0
//   Group B: {3,4,5} with cycle 3 → 4 → 5 → 3
//
// The function assumes a valid triple contains two elements from one group (dominant group)
// and one element from the other group.
//
// Behavior:
// - Validates input (all values must be in [0..5] and all must be distinct)
// - Detects the dominant group (the one appearing twice)
// - Orders the two dominant elements according to their cyclic orientation
// - Places the remaining element last
//
// Returns:
// - true  if the triple is valid and was successfully reordered
// - false if the input is invalid (out-of-range values or duplicates)
//
// The reordering is done in-place.
    // ---- validity check ----
    for (int v : t)
        if (v < 0 || v > 5)
            return false;

    if (t[0] == t[1] || t[0] == t[2] || t[1] == t[2])
        return false;

    // ---- determine dominant group ----
    bool lowGroup =
        ((t[0] < 3) + (t[1] < 3) + (t[2] < 3)) >= 2;

    int a, b, c;

    if ((t[0] < 3) == lowGroup && (t[1] < 3) == lowGroup)
    {
        a = t[0]; b = t[1]; c = t[2];
    }
    else if ((t[0] < 3) == lowGroup && (t[2] < 3) == lowGroup)
    {
        a = t[0]; b = t[2]; c = t[1];
    }
    else
    {
        a = t[1]; b = t[2]; c = t[0];
    }

    // ---- enforce cyclic orientation ----
    // maps: 0-1-2 and 3-4-5 onto same cycle via mod 3
    if (((a % 3) + 1) % 3 != (b % 3))
        std::swap(a, b);

    t = {a, b, c};
    return true;
}
Vec3D TriangularPrismBuilder::TripleNormal(triple& t) {
// Computes the normal vector of a triangular face defined by a triple,
// taking periodic boundary conditions (PBC) into account using the minimum image convention.
//
// The triangle is formed from vertices:
//   t[0] -> reference vertex
//   t[1], t[2] -> other two vertices

// Assumes:
// - Orthorhombic periodic box (diagonal box matrix)
// - m_pBox contains box lengths in each dimension
// - m_Vertices contain valid vertex pointers
//
// Note: This implementation is not valid for triclinic (tilted) simulation boxes.

    Vec3D dx1 = m_Vertices[t[1]]->GetPos() - m_Vertices[t[0]]->GetPos();
    Vec3D dx2 = m_Vertices[t[2]]->GetPos() - m_Vertices[t[0]]->GetPos();

    Vec3D box = *m_pBox;

    for (int i = 0; i < 3; i++)
    {
        if (dx1(i) >  box(i) / 2) dx1(i) -= box(i);
        if (dx1(i) < -box(i) / 2) dx1(i) += box(i);

        if (dx2(i) >  box(i) / 2) dx2(i) -= box(i);
        if (dx2(i) < -box(i) / 2) dx2(i) += box(i);
    }

    dx1.normalize();
    dx2.normalize();

    Vec3D normal = dx1 * dx2; // cross product
    normal.normalize();
    
    return normal;
}
bool TriangularPrismBuilder::UpdatePrismNormals(TriangularPrism& tp)
{
    tp.VNormals.clear();
    tp.VNormals.reserve(tp.VTriples.size());

    for (auto& tri : tp.VTriples)
    {
        tp.VNormals.push_back(TripleNormal(tri));
    }

    return true;
}
bool TriangularPrismBuilder::MakePrismMaps() {
    /**
 * @brief Builds prism topology maps either from a file or using default generation
 * 
 * This function attempts to load prism topology definitions from a user-specified
 * file. If no filename is provided, it falls back to generating default prism maps.
 * Each map in the file is read sequentially using ReadPrismMap() and stored in
 * the m_TopologyMaps container.
 * 
 * @return true  if maps were successfully loaded or default maps were created
 * @return false if the specified topology file exists but cannot be opened
 * 
 * @note The function expects the file format compatible with ReadPrismMap()
 * @see ReadPrismMap(), MakeDefaultPrismMaps()
 */
    
    
    // If a topology file is specified, attempt to load from it
    int total_map_inFile = 0;
    if (!m_TopologyFilename.empty()) {
        std::ifstream mapFile(m_TopologyFilename);
        
        if (!mapFile.is_open()) {
            std::cerr << "---> Error: (id 139565) failed to open map file: " 
                      << m_TopologyFilename << std::endl;
            return false;
        }
        
        std::string line;
        while (std::getline(mapFile, line)) {
            // Skip empty lines or comments 
            if (line.empty() || line[0] == '#') {
                continue;
            }
            
            // Read a single prism map from the file stream
            // Note: Map ID and other metadata are handled inside ReadPrismMap()
            TriangularPrism prism = ReadPrismMap(mapFile);
            total_map_inFile++;
            if(CheckPrismMapOrientation(prism)){
                m_TopologyMaps.push_back(std::move(prism));
            }
        }
        
            std::cout << "---> Note: The fusion map contains " << total_map_inFile << " fusion map";
            if (m_TopologyMaps.size() == total_map_inFile) {
                    std::cout << ", all maps will be considered.";
            } 
            else {
                std::cout << ", but only " << m_TopologyMaps.size() << " will be considered.\n";
            }
        
        return true;
    } 
    
    // No filename provided - generate default topology maps
    return MakeDefaultPrismMaps();
}
TriangularPrism TriangularPrismBuilder::ReadPrismMap(std::istream& input) {
/**
 * Reads a triangular prism definition from a stream.
 *
 * Expects exactly 6 triples of integer vertex indices.
 * Each triple is reordered via ReorderTriple() before being stored.
 *
 * @param input Input stream containing prism data.
 * @return A populated TriangularPrism.
 * @throws std::runtime_error if an invalid triple is encountered.
 */
    
    TriangularPrism prism;
    prism.VTriples.reserve(6);
    std::string clearline;

    for (int i = 0; i < 6; ++i)
    {
        int a, b, c;
        input >> a >> b >> c;
        getline(input, clearline);
        triple t = make_triple(a, b, c);

        if (!ReorderTriple(t))
        {
            throw std::runtime_error("Invalid prism triple.");
        }

        prism.VTriples.push_back(t);
    }

    return prism;
}
bool TriangularPrismBuilder::MakeDefaultPrismMaps()
{
    //==== map 1
    TriangularPrism map1;
    std::vector<triple> VTriples;

    triple t;

    t = make_triple(0, 3, 4);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    t = make_triple(0, 4, 1);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    t = make_triple(2, 1, 4);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    t = make_triple(2, 4, 5);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    t = make_triple(0, 2, 3);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    t = make_triple(2, 5, 3);
    if (!ReorderTriple(t)) return false;
    VTriples.push_back(t);

    map1.VTriples = std::move(VTriples);

    m_TopologyMaps.push_back(map1);

    return true;
}
bool TriangularPrismBuilder::CheckPrismMapOrientation(TriangularPrism & tp){
/**
 * brief: Checks whether the triangular prism face map can be made with consistent orientation.
 *
 * This function verifies that all directed edges in the prism surface occur in
 * balanced opposite pairs, ensuring consistent face winding across adjacent
 * triangles.
 *
 * For each triangle in the prism, directed edges are extracted according to the
 * stored vertex ordering. Only the edges (v2 → v3) and (v3 → v1) are considered,
 * since the (v1 → v2) edge is shared with the original mesh and handled elsewhere.
 *
 * The method counts occurrences of each directed edge and confirms that every
 * edge has a corresponding opposite-directed edge with the same multiplicity.
 * If any imbalance is found, the prism orientation is considered invalid.
 *
 * @param tp Triangular prism structure containing vertex triples (VTriples).
 * @return true if all edges are consistently oriented; false otherwise.
 */
    std::vector<triple> VTriples = tp.VTriples;
    
     // Count occurrences of each directed edge
    std::unordered_map<PairMap, int, PairMapHash> edgeCount;
    
    // Count directed edges (only from vertex 2->3 and 3->1 as per original)
    for (auto tmap : VTriples) {
        edgeCount[{tmap[1], tmap[2]}]++;  // edge from vertex 2 to 3
        edgeCount[{tmap[2], tmap[0]}]++;  // edge from vertex 3 to 1
    }
    
    // Verify each directed edge has its opposite (C++14 compatible)
    for (const auto& entry : edgeCount) {
        const PairMap& edge = entry.first;
        int count = entry.second;
        
        PairMap opposite = {edge[1], edge[0]};
        auto it = edgeCount.find(opposite);
        
        // Must exist and have same count (balanced pairs)
        if (it == edgeCount.end() || it->second != count) {
            return false;
        }
    }
    
    return true;    
}