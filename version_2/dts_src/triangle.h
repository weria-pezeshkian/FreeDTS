#if !defined(AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_)
#define AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "Vec3D.h"

/*
 * @file triangle.h
 * @brief Definition of the Triangle class, representing a geometric triangle in 3D space.
 * This class encapsulates functionality related to triangles, such as storing vertex pointers,
 * calculating area and normal vectors, and updating triangle properties.
 * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
 */
class vertex;
class triangle {
public:
    triangle(int id);
    triangle(int id, vertex *v1, vertex *v2, vertex *v3);
    ~triangle();

    // Getters
    inline const int GetTriID() const { return m_ID; }
    inline vertex *GetV1() { return m_V1; }
    inline vertex *GetV2() { return m_V2; }
    inline vertex *GetV3() { return m_V3; }
    inline double GetArea() { return m_Area; }
    inline double GetVolume() { return m_Volume; }
    inline Vec3D GetAreaVector() { return m_AreaVector; }
    inline Vec3D GetNormalVector() { return m_Normal; }

    // Only needed for visualization
    inline bool GetRepresentation() { return m_Representation; }

public:
    // Update functions
    void UpdateRepresentation(bool); // For visualization output
    void UpdateNormal_Area(Vec3D *Box); // May be sent to other classes
    Vec3D CalculateNormal(Vec3D Box); // important note: This function does not change the NormalVector member variable
    void UpdateNormal_Area(Vec3D& norm, double& area); // sets normal and area
    void UpdateVertex(vertex *v1, vertex *v2, vertex *v3); // Update vertices
    void UpdateVolume(double vol);
    void UpdateID(int id); // Should not be used for active triangles

    // Copy and reverse state
    void Copy();
    void Reverse2PreviousCopy();
    void ConstantMesh_Copy();
    void ReverseConstantMesh_Copy();
private:
    // Original state members
    int m_ID;
    bool m_Representation;
    vertex *m_V1;
    vertex *m_V2;
    vertex *m_V3;
    Vec3D m_Normal;
    Vec3D m_AreaVector;
    double m_Area;
    double m_Volume;

    // Old state members to maintain previous state
    vertex *m_oldV1;
    vertex *m_oldV2;
    vertex *m_oldV3;
    Vec3D m_oldNormal;
    Vec3D m_oldAreaVector;
    double m_oldArea;
    double m_oldVolume;
    
private:
    double adjust_periodic(double d, double box_dim);
};

#endif
