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
struct pot_triangle {    // data structure for a potential triangle
        int id;
        int cid; // connected id
        vertex *pv1;
        vertex *pv2;
        vertex *pv3;
        links *pl1;
        links *pl2;
        links *pl3;
    };
    struct pair_pot_triangle {    // data structure for a potential triangle
        bool state;
        int id;
        pot_triangle PT1;
        pot_triangle PT2;
        std::vector <links *> ConnectingLinks;       // this does not include the mirror links
        std::vector <triangle *> ConnectingTriangles;

    };
class Three_Edge_Scission :  public AbstractDynamicTopology,
                             public AbstractSimulation,
                             public MESH { // to use for polymorphism

public:
    Three_Edge_Scission(int period, State *pState);
    ~Three_Edge_Scission();
    void Initialize();
    bool MCMove(int step);
    std::string CurrentState();

    inline  std::string GetDerivedDefaultReadName() {return "Three_Edge_Scission";}
    inline static std::string GetDefaultReadName() {return "Three_Edge_Scission";}
    
private:
    bool MCScissionMove(int step);
    bool MCFussionMove(int step);
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    template<typename T>
    bool AddtoVectorCarefully(T* z, std::vector<T*> &vect);
    
    bool Anglevalid4Vhole(vertex *v1, double minangle);
    bool CorrectOrientation(pot_triangle &p1,pot_triangle &p2);
    pair_pot_triangle connected_2pot_triangles(pot_triangle potT1, pot_triangle potT2);
    std::vector<pair_pot_triangle> FindPotentialTriangles();
    std::vector <triangle *> DoAScission(pair_pot_triangle &pair);
    bool ReverseAScission(pair_pot_triangle &pair, triangle *t1, triangle *t2);   // this is the exact reverse action of DoAScission; different from DoAFussion

    bool DoAFussion(pair_pot_triangle pair);
    triangle * CreateATriangleFromAPotentialTriangle(pot_triangle &p1);
    
                                 
    int m_Period;
    State *m_pState;




template<typename T>
void KeepOneOccurrence(std::vector<T*> &vect);


private:
bool do_Simulation(){
    std::cout<<" ---> error, 0009989 this should have been called \n";
    return false;
}
    
};
#endif
