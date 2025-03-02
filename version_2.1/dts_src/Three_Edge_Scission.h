#if !defined(THREE_EDGE_SCISSION_H_INCLUDED)
#define THREE_EDGE_SCISSION_H_INCLUDED

#include "AbstractDynamicTopology.h"
#include "AbstractSimulation.h"
#include "SimDef.h"
#include "Vec3D.h"
#include "MESH.h"

class State;
class triangle;
class vertex;
class links;

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
    int cid;      ///< ID of another potential triangle that is connected to this one, if any (-1 if not connected).

};
struct pair_pot_triangle {    // data structure for a pair of potential_triangle
        bool state;
        int id;
        pot_triangle PT1;
        pot_triangle PT2;
        std::vector <links *> ConnectingLinks;       // this does not include the mirror links
        std::vector <triangle *> ConnectingTriangles;

};
struct fussion_site {    // data structure for a pair of triangle
    links* l1;
    links* l2;
    double dist[3][3];
    int no_conf;
        
};
class Three_Edge_Scission :  public AbstractDynamicTopology { // to use for polymorphism

public:
    Three_Edge_Scission(int period, State *pState);
    ~Three_Edge_Scission();
    void Initialize();
    bool MCMove(int step);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "Three_Edge_Scission";}
    inline static std::string GetDefaultReadName() {return "Three_Edge_Scission";}

    
    
    
private:
    //--- functions for fission
    bool CorrectOrientation(pot_triangle &p1,pot_triangle &p2);
    bool CheckFaceAngle(pair_pot_triangle &pair);
    bool Is_A_Neck(pot_triangle potT1, pot_triangle potT2, pair_pot_triangle & neck);
    std::vector<pair_pot_triangle> FindNecks();
    std::vector <triangle *> DoAScission(pair_pot_triangle &pair);
    bool ReverseAScission(pair_pot_triangle &pair, triangle *t1, triangle *t2);   // this is the exact reverse action of DoAScission; different from DoAFussion
    std::vector<links*> GetEdgesWithInteractionChange(pair_pot_triangle &pair);
    triangle * CreateATriangleFromAPotentialTriangle(pot_triangle &p1);
    bool ScissionByMC(pair_pot_triangle &pair_t, double thermal);
    
    //--- functions for fussions
    bool FussionByMove(fussion_site &pair_tri, double thermal);
    bool CheapScane(links *l1, links *l2, fussion_site &p_T);
    bool BuildScane(fussion_site &p_T);

    std::vector<fussion_site> FindPotentialFussionSites();


private:


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
    void ActivateAGhostTriangle(triangle *tr);

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
    
};
#endif
