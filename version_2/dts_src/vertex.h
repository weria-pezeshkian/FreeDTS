#if !defined(AFX_vertex_H_9P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_vertex_H_9P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_
#ifdef _OPENMP
# include <omp.h>
#endif
#include "Voxel.h"
#include "SimDef.h"
#include "Vec3D.h"
#include "Tensor2.h"
#include "inclusion.h"
#include "VertexVectorFields.h"
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Vertex object: can be defined by and integer and three double numbers for the positions of the vertices
 Also there are a few lists pointing to the nighbouring links and trinagles and also the inclusion if the vertex own one.
 *******************/
class links;
class triangle;
class MESH;
class vertex : public VertexVectorFields {
    
public:
	vertex(MESH* pMesh, int id, double x, double y, double z);
	vertex(MESH* pMesh, int id);
    vertex();

	 ~vertex();

// A set of function to get vertex variables (since all of them are private)
	    inline int GetVID()                                 {return m_ID;}
        inline int GetVertexType()                          {return m_VertexType;}
        inline int  GetDomainID()                           {return m_DomainID;}
        inline Tensor2  GetL2GTransferMatrix()              {return m_T_Local_2_Global;}
        inline Tensor2  GetG2LTransferMatrix()              {return m_T_Global_2_Local;}
        inline inclusion* GetInclusion()                    {return m_pInclusion;}
        inline bool VertexOwnInclusion()                    {return m_OwnInclusion;}
        inline Vec3D GetNormalVector()                      {return m_Normal;}
        inline double GetEnergy()                           {return m_Energy;}
        inline Vec3D *GetBox()                              {return m_pBox;}
        inline int GetGroup()                               {return m_Group;}
        inline std::string GetGroupName()                   {return m_GroupName;}
        inline std::vector <links *> GetVLinkList()             {return m_VLinkList;}
        inline std::vector <triangle *> GetVTraingleList()         {return m_VTraingleList;}
        inline std::vector <vertex *> GetVNeighbourVertex()     {return m_VNeighbourVertex;}
        inline Voxel<vertex> * GetVoxel()     {return m_pVoxel;}

//---> functions to access positions
        inline double GetVXPos()                            {return m_X;}
        inline double GetVYPos()                            {return m_Y;}
        inline double GetVZPos()                            {return m_Z;}
        inline double GetXPos()                             {return m_X;}
        inline double GetYPos()                             {return m_Y;}
        inline double GetZPos()                             {return m_Z;}
        inline  Vec3D GetPos()                              {return   Vec3D(m_X,m_Y,m_Z);}

//---->
        //---> for if the vertex is surface vertex
        inline double GetP1Curvature()              {return m_PrincipalCurvature_1;}// surface curvature
        inline double GetP2Curvature()              {return m_PrincipalCurvature_2;}// surface curvature
        inline double GetArea()                     {return m_Area;}
        //---> for if the vertex is edge vertex
        inline double GetGeodesicCurvature()        {return m_Geodesic_Curvature;}// surface curvature
        inline double GetNormalCurvature()          {return m_Normal_Curvature;}// surface curvature
        inline double GetLength()                   {return m_VLength;}// surface curvature
        inline links* GetPrecedingEdgeLink()        {return m_pPrecedingEdgeLink;}// preceding link at the edge
        inline links* GetEdgeLink()                 {return m_pEdgeLink;}// preceding link at the edge

    

public:
    
    bool SetCopy();            // Copies the key ellements into the old type
    bool Reverse2PreviousCopy();  // reverse the edge to the value set at the time of MakeCopy()
    
    void ConstantMesh_Copy();
    void ReverseConstantMesh_Copy();
    
    void EnergyCopy();
    void ReverseEnergyCopy();
    
  // A set of functions to update vertex variables
  void UpdateVXPos(double x);    // a function for update the x position of a vertex
  void ScalePos(double lx, double ly, double lz);    // a function for update the x position of a vertex
  void UpdateVYPos(double y);   // a function for update the y position of a vertex
  void UpdateVZPos(double z);   // a function for updates the z position of a vertex
  void PositionPlus(double dx, double dy, double dz);   // a function to increase the positions by dx, dt, dz
  void UpdateP1Curvature(double p1_curvature);
  void UpdateP2Curvature(double p2_curvature);
  void UpdateGroupName(std::string z); // A vertex can only have one group name, is different from group id
  void UpdateBox(Vec3D *z);
  void UpdateEnergy(double);
  void UpdateNormal_Area(Vec3D,double); // vis
  void UpdateL2GTransferMatrix(Tensor2 v);
  void UpdateG2LTransferMatrix(Tensor2 v);
  void UpdateInclusion(inclusion* );
  void UpdateOwnInclusion(bool );
  void UpdateVID(int); // this should never be called, only for temporarily vis functions
  void AddtoLinkList(links* z);
  bool AddtoLinkListCarefully(links* z);
  void AddtoTraingleList(triangle * z);
  bool AddtoTriangleListCarefully(triangle* t);
  void AddtoNeighbourVertex(vertex* z);
  bool AddtoNeighbourVertexCarefully(vertex* v);
  void RemoveFromLinkList(links* z);
  bool RemoveFromLinkListCarefully(links* l);
  void RemoveFromTraingleList(triangle * z);
  void RemoveFromNeighbourVertex(vertex* z);
  void UpdateGroup(int z);
  void UpdateVoxel(Voxel<vertex> * pVoxel);
  void UpdateDomainID(int domain_id);

public:
    bool CheckVoxel();
    bool UpdateVoxelAfterAVertexMove(); // update the voxel after the move has happened.
    double SquareDistanceFromAVertex(vertex* pv2);
    double SquareDistanceOfAVertexFromAPoint(double X, double Y, double Z, vertex* pv2);
    bool UpdateVFGlobalDirectionFromLocalDirection();   // check if the global direction of the vector fields can be obtianed from local
    bool UpdateVFLocalDirectionFromGlobalDirection();   // check if the global direction of the vector fields can be obtianed from local
    bool ReverseVFLocalDirection();

    
    friend std::ostream& operator<<(std::ostream& os, const vertex& obj);
    

private:

    int m_ID;         // ID of the vertex, a unique number
    int m_DomainID;
    double m_X;       // X coordinate
    double m_Y;       // Y coordinate
    double m_Z;       // Z coordinate
    std::vector <triangle *> m_VTraingleList;  // A list holding all the nighbouring triangles
    std::vector <links *> m_VLinkList;      // A list holding all the nighbouring edges (linkss)
    std::vector <vertex *> m_VNeighbourVertex; // A list holding all the nighbouring vertexes
    inclusion *m_pInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    bool m_OwnInclusion;                        // to check if the vertex own any inclusion
    double m_Area;                              // area of the vertex
    int m_Group;            // Id of a group that the vertex belong too
    private:
    Vec3D m_Normal;
    Vec3D *m_pBox;
    double m_Energy;
    Tensor2  m_T_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_T_Global_2_Local;        //  global to local transformation matrix
    std::string m_GroupName;
     Voxel<vertex> * m_pVoxel;
    double m_PrincipalCurvature_1;
    double m_PrincipalCurvature_2;
    double m_MeanCurvature;
    double m_GaussianCurvature;
    
public:
    // lets have them public for now
    //================================================
    // development of Aug 2023; for membranes with hole
    //a*kg^2+b*k^2    note: any spontaneous one will come from inclusions
    
    // we might need to define n and t and b as well, for now lets ignore it
    //================================================
    double m_Geodesic_Curvature;          // Edge Vertex Curvature
    double m_Normal_Curvature;          // Edge Vertex Curvature
    double m_VLength;                       // length of the vertex
    
    int m_VertexType;                   // 0 surface vertex; 1 edge vertex;
    links * m_pEdgeLink;
    links * m_pPrecedingEdgeLink;// preceding link at the edge

    
// members for copying.
private:
    double m_OldX;       // X coordinate
    double m_OldY;       // Y coordinate
    double m_OldZ;       // Z coordinate
    std::vector <triangle *> m_OldVTraingleList;  // A list holding all the nighbouring triangles
    std::vector <links *> m_OldVLinkList;      // A list holding all the nighbouring edges (linkss)
    std::vector <vertex *> m_OldVNeighbourVertex; // A list holding all the nighbouring vertexes
    inclusion *m_OldpInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    bool m_OldOwnInclusion;                        // to check if the vertex own any inclusion
    double m_OldArea;                              // area of the vertex
    int m_OldGroup;            // Id of a group that the vertex belong too
    Vec3D m_OldNormal;
    double m_OldEnergy;
    Tensor2  m_OldT_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_OldT_Global_2_Local;        //  global to local transformation matrix
    std::string m_OldGroupName;
    Voxel<vertex> * m_OldpVoxel;
    double m_OldGeodesic_Curvature;          // Edge Vertex Curvature
    double m_OldNormal_Curvature;          // Edge Vertex Curvature
    double m_OldVLength;                       // length of the vertex
    int m_OldVertexType;                   // 0 surface vertex; 1 edge vertex;
    links * m_OldpEdgeLink;
    links * m_OldpPrecedingEdgeLink;// preceding link at the edge
    MESH* m_pMesh;
    double m_OldPrincipalCurvature_1;
    double m_OldPrincipalCurvature_2;
    double m_OldMeanCurvature;
    double m_OldGaussianCurvature;
    
// OpenMP methods Oct 2024    
#ifdef _OPENMP
public:
bool CheckLockVertex();
void LockVertex();
bool UnlockVertex();
void LockNeighbourVertex();
void UnlockNeighbourVertex();
bool CheckLockNeighbourVertex();
    static bool CheckLockVectorVertex(const std::vector<vertex*>&);
    static void UnockVectorVertex(std::vector <vertex *>);
private:
omp_lock_t m_Lock;
#endif
};


#endif
