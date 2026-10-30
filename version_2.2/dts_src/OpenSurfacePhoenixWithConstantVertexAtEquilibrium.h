#if !defined(AFX_OpenSurfacePhoenixWithConstantVertexAtEquilibrium_H_INCLUDED_)
#define AFX_OpenSurfacePhoenixWithConstantVertexAtEquilibrium_H_INCLUDED_
#include "AbstractOpenEdgeEvolution.h"
#include "AbstractSimulation.h"
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 
 This is a generic command to create open surfaces. 
 
 The four main the the most important functions are
 
    1) links* CreateALink(vertex *);
    2) bool KillALink(links *);
    3) triangle* CloseATriangleHole(vertex *v1);
    4) bool KillATriangle(links *l1);
 
 
 */

class State;
class OpenSurfacePhoenixWithConstantVertexAtEquilibrium : public AbstractOpenEdgeEvolution { // to use for polymorphism
public:
    OpenSurfacePhoenixWithConstantVertexAtEquilibrium(int period, double rate, State *pState);
    ~OpenSurfacePhoenixWithConstantVertexAtEquilibrium();

public:
    void Initialize();
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "OpenSurfacePhoenix";}
    inline static std::string GetDefaultReadName() {return "OpenSurfacePhoenix";}
    
    bool Move(int step); // this is the main function to perform the move

    
private:
    bool MCAttemptedToAddALink();
    bool MCAttemptedToRemoveALink();


    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoVertexList(vertex* z, std::vector<vertex*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    bool Linkisvalid(vertex *);



    // the main hard part of the code. 
private:  // this functions could be in princeple public, but no need for now
    
    links* CreateALink(vertex *);    // creates a link between two connected edge at the edge
    bool KillALink(links *);         // kills an edge that is on the edge. not an edge on the surface
    triangle* CloseATriangleHole(vertex *v1);  //if an edge contains only three links, it will close it
    bool KillATriangle(links *l1);      // kills a triangle anywhere

    
    
    // helper functions
    int GetSurfVFirstRingSize() const; // Computes the number of unique vertices in the first ring of the surface.
    
private:
    State *m_pState;
    int m_Period;
    double m_NumberOfMovePerStep;   // how many updates should be made per step
    
    std::vector<links*>&          m_pEdgeL;
    std::vector<links*>&          m_pGhostL;
    std::vector<links*>&          m_pRightL;
    std::vector<links*>&          m_pLeftL;
    std::vector<links*>&          m_pActiveL;
    std::vector<vertex*>&         m_pEdgeV;
    std::vector<vertex*>&         m_pSurfV; // all the active vertices  surf
    std::vector<triangle*>&       m_pGhostT;
    std::vector<triangle*>&       m_pActiveT;
    MESH *m_pMesh;
    Vec3D *m_pBox;
    
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;
    int &m_No_VectorFields_Per_V;
private:
    bool do_Simulation(){
        std::cout<<" ---> error, 999o1o this should have been called \n";
        return false;
    }
};
#endif
