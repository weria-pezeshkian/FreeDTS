/**
 * @class TopologyChangeByTriangularPrism
 * @brief Monte Carlo dynamic-topology algorithm based on triangular-prism
 *        topology transformations.
 *      Started June 2024
 *      End: June 2026
 *
 *   Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 *   Copyright (c) Weria Pezeshkian
 *
 * This class implements topology-changing moves for triangulated membrane
 * surfaces by representing local fusion and scission events as transformations
 * of a six-vertex triangular prism. The method allows the mesh connectivity
 * and surface genus to evolve while preserving a valid triangulation.
 *
 * The algorithm supports:
 * - Detection of neck-like structures that can undergo scission.
 * - Detection of nearby triangle pairs that can undergo fusion.
 * - Construction of admissible prism topologies using a
 *   TriangularPrismBuilder.
 * - Execution and reversal of local topology modifications.
 * - Monte Carlo acceptance/rejection using the Metropolis criterion.
 * - Consistent updates of mesh geometry, curvature, inclusion interactions,
 *   vector-field interactions, and global constraint energies.
 *
 * Topology changes are performed using preallocated ghost triangles and
 * ghost links, avoiding dynamic memory allocation during simulation.
 * Rejected moves are exactly reversible through stored local copies of
 * affected vertices, links, and interaction energies.
 *
 * Energy contributions considered during acceptance may include:
 * - Local vertex energies
 * - Inclusion interaction energies
 * - Vector-field interaction energies
 * - Volume constraint energy
 * - Total area constraint energy
 * - Global curvature energy
 *
 * The class operates periodically during the simulation and updates the
 * surface genus according to the current Euler characteristic of the mesh.
 *
 * Assumptions:
 * - The mesh is represented as an oriented triangular surface.
 * - Sufficient ghost triangles and ghost links are available for topology
 *   modifications.
 * - Local topology transformations preserve mesh validity and orientation.
 *
 *
 * @see TriangularPrismBuilder
 * @see AbstractDynamicTopology
 */
#if !defined(TOPOLOGYCHANGEBYTRIANGULARPRISM_H_INCLUDED)
#define TOPOLOGYCHANGEBYTRIANGULARPRISM_H_INCLUDED
#include <array>
#include <unordered_set>
#include "AbstractDynamicTopology.h"
#include "AbstractSimulation.h"
#include "SimDef.h"
#include "Vec3D.h"
#include "MESH.h"
#include "Tensor2.h"
#include "TriangularPrismBuilder.h"

class State;
class triangle;
class vertex;
class links;


struct PotTriangleKey {
    vertex* v1;
    vertex* v2;
    vertex* v3;

    bool operator==(const PotTriangleKey& other) const {
        return v1 == other.v1 && v2 == other.v2 && v3 == other.v3;
    }
};

struct PotTriangleKeyHash {
    size_t operator()(const PotTriangleKey& k) const {
        size_t h1 = std::hash<vertex*>()(k.v1);
        size_t h2 = std::hash<vertex*>()(k.v2);
        size_t h3 = std::hash<vertex*>()(k.v3);

        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};
struct pot_triangle {
    /**
     * @brief Data structure representing a potential triangle in the mesh.
     * This structure represents a potential triangle in the mesh. It consists of three vertices and their corresponding edges.
     * Note that these edges do not yet form a any triangle in the mesh.
     */
    int id;       ///< Unique identifier for the potential triangle.
    vertex *pv1;  ///< Pointer to the first vertex of the potential triangle.
    vertex *pv2;  ///< Pointer to the second vertex of the potential triangle.
    vertex *pv3;  ///< Pointer to the third vertex of the potential triangle.
    links *pl1;   ///< Pointer to the link (edge) between pv1 and pv2.
    links *pl2;   ///< Pointer to the link (edge) between pv2 and pv3.
    links *pl3;   ///< Pointer to the link (edge) between pv3 and pv1.

};
struct fission_site {    // data structure for a pair of potential_triangle
    std::array<vertex*, 3> v_ver;
    std::array<vertex*, 3> u_ver;
    std::vector<links*> C_Links;   // connecting links, this does not include the mirror links
    std::array<links*, 3> V_links; // Top 3 links, this does not include the mirror links
    std::array<links*, 3> U_links; // Top 3 links, this does not include the mirror links
};
struct PairHash {
    size_t operator()(const std::pair<triangle*, triangle*>& p) const {
        size_t h1 = std::hash<triangle*>()(p.first);
        size_t h2 = std::hash<triangle*>()(p.second);
        return h1 ^ (h2 << 1);
    }
};
struct fusion_site {    // data structure for a fusion site
    triangle* t1;
    triangle* t2;
    TriangularPrism topology;     
};
struct fusion_outcome {    // data structure for a fusion site
    std::vector<triangle*> newTriangles;
    std::vector<links*> newLinks;
    std::vector<links*> oldLinks;
    std::vector<links*> pHalfnewLinks;
    std::vector<vertex*> vertices;
    triangle* t1;
    triangle* t2;
};
class TopologyChangeByTriangularPrism :  public AbstractDynamicTopology { // to use for polymorphism

public:
    TopologyChangeByTriangularPrism(std::string inputdata, State *pState);
    ~TopologyChangeByTriangularPrism();
    void Initialize();
    bool MCMove(int step);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "TriangularPrism";}
    inline static std::string GetDefaultReadName() {return "TriangularPrism";}

    
    
    
private:
    //--- functions for fission
    bool CheckFaceAngleForFissionSite(fission_site &f_site);
    //std::vector<pair_pot_triangle> FindNecks();
    std::vector <triangle *> DoAScission(fission_site &pair);
    bool ReverseAScission(fission_site &pair);   // this is the exact reverse action of DoAScission; different from DoAFussion
    std::vector<links*> GetEdgesWithInteractionChange(fission_site &f_site);
    triangle * CreateATriangleFromAPotentialTriangle(pot_triangle &p1);
    bool ScissionByMC(fission_site &pair_t, double thermal,int n_fi_site);
    
    //--- functions for fussions
    bool FusionByMove(fusion_site &pair_tri, double thermal, int no_cur_fu_sites);


    std::vector<fusion_site> FindPotentialFusionSites();

// fusion functions
    bool FusionSites_AreNotNeighbours(triangle *t1, triangle *t2);
    bool Fuse_MeshViaTwoTriangles(fusion_site &pair_tri, fusion_outcome &data);
    bool Reverse_Fuse_MeshViaTwoTriangles(fusion_outcome &fusion_mesh);
    std::vector<links*> Get_EdgesFusionAffect(std::vector<vertex*> &Vver);
// fission functions
    bool Is_A_Neck(pot_triangle potT1, pot_triangle potT2, fission_site & neck);
    std::vector<fission_site> FindNecks();
    triangle * CreateTriangleByTriple(std::array<vertex*, 3>& v_ver, std::array<links*, 3>& v_links);
private:
   bool  VoxelizeTriangles(double voxsize);
  // bool (*AcceptanceProbability)(double dE);

private:
template<typename T>
void KeepOneOccurrence(std::vector<T*> &vect);
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    template<typename T>
    bool AddtoVectorCarefully(T* z, std::vector<T*> &vect);
    
    
    void MakeALinkGhost(links *l);
    void ActivateAGhostLink(links *l);
    void MakeATriangleGhost(triangle *tr);

private:    
    std::vector<links*>&          m_pEdgeL;
    std::vector<links*>&          m_pGhostL;
    std::vector<links*>&          m_pRightL;
    std::vector<links*>&          m_pLeftL;
    std::vector<links*>&          m_pActiveL;
    std::vector<vertex*>&         m_pEdgeV;
    std::vector<vertex*>&         m_pSurfV; // all the active vertices  surf
    std::vector<triangle*>&       m_pGhostT;
    std::vector<triangle*>&       m_pActiveT;
    Vec3D &m_Box;
    
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;
    int &m_No_VectorFields_Per_V;
    int m_Period;
    State *m_pState;
    Voxelization<triangle>  *m_pTriVoxelization;
    TriangularPrismBuilder *m_pTriangularPrismBuilder;
    
    bool m_No_Vectorfield;
    double m_Local_KG; // we have this because fission and fusion has discrete contribution to gaussain term, so, we do not need to calcute the energy for it. 
private:
std::string m_StreamInputs;    
std::string m_PrismMapTopologyFile;    

};
#endif
