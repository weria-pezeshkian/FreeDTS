#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include <map>
#include "SimDef.h"
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "CreateMashBluePrint.h"
/*
 * @class MESH
 * @brief This class provides functionalities for manipulating the mesh structure.
 *
 * The MESH class contains methods for creating, updating, and manipulating the mesh.
 * It handles various mesh elements including vertices, triangles, and links. It also
 * provides methods for handling inclusions, groups, and boundary conditions.
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */
class MESH
{
public:
    MESH();  // Default constructor for MESH.
    ~MESH();  // Destructor for MESH.


    inline std::vector<vertex*>&        GetActiveV()               {return m_pActiveV;}
    inline std::vector<triangle*>&      GetActiveT()               {return m_pActiveT;}
    inline std::vector<links*>&      GetActiveL()               {return m_pActiveL;}
    inline std::vector<vertex*>&           GetSurfV()                 {return m_pSurfV;}
    inline std::vector<vertex*>&           GetEdgeV()                 {return m_pEdgeV;}
    inline std::vector<links*>&      GetEdgeL()                      {return m_pEdgeL;}
    inline std::vector<links*>&      GetRightL()                    {return m_pHL;}
    inline std::vector<links*>&        GetLeftL()                 {return m_pMHL;}
    inline std::vector<links*>&        GetGhostL()                 {return m_pGhostL;}
    inline std::vector<triangle*>&     GetGhostT()                {return m_pGhostT;}

    inline std::vector<inclusion*>&  GetInclusion()             {return m_pInclusion;}
    std::map<std::string, std::vector<vertex*> > &GetGroups()    {return m_Groups;}


    inline Vec3D                   *GetBox()                        {return m_pBox;}
    inline Vec3D&                   Link2ReferenceBox()             {return m_Box;}
    inline const bool               GetHasCrossedPBC()  const       {return m_MeshCrossedPBC;}
    std::vector <InclusionType*>    GetInclusionType()     const    {return m_pInclusionType;}
    inline  int&                    GetNoVFPerVertex()         {return m_No_VectorFields_Per_V;} // returning the number of vector field per vertex

    
    
    
    inline void UpdateCrossedPBC(bool newValue){
        if(!m_MeshCrossedPBC)
            m_MeshCrossedPBC = newValue;
        return;
    }
    
// some helpful static functions
    static double SquareDistanceBetweenTwoVertices(vertex *p_v1, vertex* p_v2, Vec3D Box);
    double SquareDistanceBetweenTwoVertices(vertex * v1,vertex * v2);



public:
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    
// -- it could be optimised
    void MakeALinkGhost(links* l);
    void MakeATriangleGhost(triangle* tri);
    // this has not been completed yet
    bool Flip(links* l);
    bool MakeAVertexGhost(vertex* v);
    bool EdgeV2SurfV(vertex* v);
//----
    bool UpdateGroupFromIndexFile(std::string &indexfilename);
    void UpdateNoVFPerVertex(int number);
    void CenterMesh();    // this function centers the mesh inside the box. For broken systems it may not work nicely
    bool GenerateMesh(MeshBluePrint meshblueprint);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);
    std::vector <InclusionType*> m_pInclusionType;
    std::vector <InclusionType> m_InclusionType;
    bool CheckMesh(double min_l, double max_l, double min_angle, Voxelization<vertex>  *pVoxelization);

protected:
    std::vector<vertex*>        m_pActiveV; // all the active vertices edge + surf
    std::vector<vertex*>        m_pSurfV; // all the active vertices  surf
    std::vector<vertex*>        m_pEdgeV;  // edge
    std::vector<links*>         m_pActiveL;   // all the links
    std::vector<links*>         m_pHL;
    std::vector<links*>         m_pMHL;
    std::vector<links*>         m_pEdgeL;
    std::vector<triangle*>      m_pActiveT;
    std::vector<inclusion*>     m_pInclusion;
    std::vector<triangle*>      m_pGhostT; // Some trinagles for initial storing
    std::vector<links*>         m_pGhostL;
    std::vector<vertex*>        m_pGhostV;
    Vec3D                       *m_pBox;
    std::map<std::string, std::vector<vertex*> > m_Groups;


    
private:
    int m_No_VectorFields_Per_V;
    bool m_MeshCrossedPBC;
//========== this variables should be fully hidden from anything =======================
    std::vector<vertex>         m_Vertex;                           //                ||
    std::vector<triangle>       m_Triangle;                           //              ||
    std::vector<links>          m_Links;                           //                 ||
    std::vector<inclusion>      m_Inclusion;                           //             ||
    Vec3D                       m_Box;                           //                   ||
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing       ||
    std::vector<links>          m_GhostL;                           //                ||
    std::vector<vertex>         m_GhostV;                           //                ||
//======================================================================================

};



#endif
