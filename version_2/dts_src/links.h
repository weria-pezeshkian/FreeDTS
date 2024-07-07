#if !defined(AFX_links_H_6R4B21B8_K13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_links_H_6R4B21B8_K13C_5648_BF23_124095086234__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "Vec3D.h"
#include "Tensor2.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
link object (edge):
 */
class vertex;
class links
{
public:
    
	links(int id, vertex *v1, vertex *v2, triangle *t1);
	links(int id);
	 ~links();

// a set of inline function to access the private meembers
        inline const int GetID()                 const      {return m_ID;}
        inline vertex *GetV1()                   const      {return m_V1;}
        inline vertex *GetV2()                   const      {return m_V2;}
        inline vertex *GetV3()                   const      {return m_V3;}
        inline  bool GetVisualize ()             const      {return m_Show;}
        inline bool GetMirrorFlag ()             const      {return m_mirorflag;}
        inline triangle *GetTriangle()           const      {return m_T1;}
        inline links   *GetMirrorLink()          const      {return m_mirorlink;}
        inline links   *GetNeighborLink1()       const      {return m_neighborlink1;}
        inline links   *GetNeighborLink2()       const      {return m_neighborlink2;}
        inline Vec3D GetNormal()                 const      {return m_Normal;}
    	inline Vec3D GetBe()                     const     	{return m_Be;}
    	inline double GetHe()                    const      {return m_He;}
        inline double GetIntEnergy()             const      {return m_IntEnergy;}
        inline double GetEdgeSize()              const      {return m_EdgeSize;}
        inline Vec3D GetEdgeVector()             const      {return m_EdgeVector;}
        inline int GetEdgeType()                 const      {return m_LinkType;}
public:
    void UpdateTriangle(triangle *v);
    void UpdateV(vertex *v1,vertex *v2,vertex *v3);
    void UpdateV3(vertex *v3);
    void UpdateMirrorLink(links *v);
    void UpdateNeighborLink1(links *v);
    void UpdateNeighborLink2(links *v);
    void UpdateVisualize(bool v);
    void UpdateNormal();
    void PutNormal(Vec3D);
    void UpdateShapeOperator(Vec3D *);
    void PutShapeOperator(Vec3D Be,double He);
    void UpdateLinkSide(double l_side);
    void UpdateEdgeVector(Vec3D e_vector, double size);
    void Flip();
    bool Flip(links *pedge);
    bool Reverse_Flip(links *pedge);
    double GetVFIntEnergy();
    double GetVFIntEnergy(int layer);

    void UpdateMirrorFlag(bool v);
    void UpdateSimTimeStep(int v);
    void UpdateIntEnergy(double en);
    bool UpdateVFIntEnergy(int vf_id, double en);  // vector field
    void UpdateVFIntEnergy(std::vector<double> VF_EN);
    void InitializeVFIntEnergy(int no_vf);   //     //-- in each run, this function should run only once

    void UpdateEdgeVector(Vec3D *pbox);
    void PutEdgeVector(Vec3D , double);

public:
    bool CheckFaceAngleWithMirrorFace(double &minangle);  // calaculates the dot(n1,n2) of this edge trinagle with the
    bool CheckFaceAngleWithNextEdgeFace(double &minangle); // the face of the mirror
private:
    int m_ID;
    triangle *m_T1;
    vertex  *m_V1;
    vertex  *m_V2;
    vertex  *m_V3;
    links   *m_mirorlink;
    links   *m_neighborlink1; /// triangle is 1-2-3;  the link is 1->2   this point to 2->3
    links   *m_neighborlink2; /// the link is 1->2   this is 3->1
    bool m_Show;
    bool m_mirorflag;
    Vec3D m_Normal;    // average of the two trinagule normal for edge links, it is
    Vec3D m_Be;
    double m_He;
    double m_IntEnergy;
    std::vector<double> m_VFieldIntEnergy;    // interaction energy of the vector fields

    
public:
    Vec3D m_EdgeVector;    // a vector along the edge
    double m_EdgeSize;    // size of the edge
    int m_LinkType;  // 0 surface link; 1 edge link;

public:
    bool SetCopy();            // Copies the key ellements into the old type
    bool Reverse2PreviousCopy();  // reverse the edge to the value set at the time of MakeCopy()
    bool Copy_InteractionEnergy();            // Copies the the inclsuion interacion energy
    bool Reverse_InteractionEnergy();            // Reverse the the inclsuion interacion energy
    bool Copy_VFInteractionEnergy();            // Copies the vector field interactions
    bool Reverse_VFInteractionEnergy();            // Reverse the vector field interactions

    void ConstantMesh_Copy();
    void ReverseConstantMesh_Copy();
    
    // this are old value of the key variables we need this for the copying
private:
    triangle *m_OldT1;     //
    vertex  *m_OldV1;
    vertex  *m_OldV2;
    vertex  *m_OldV3;
    links   *m_Oldmirorlink;
    links   *m_Oldneighborlink1; /// triangle is 1-2-3;  the link is 1->2   this point to 2->3
    links   *m_Oldneighborlink2; /// the link is 1->2   this is 3->1
    bool m_Oldmirorflag;
    int m_OldLinkSide;
    Vec3D m_OldNormal;    // average of the two trinagule normal for edge links, it is
    Vec3D m_OldBe;
    double m_OldHe;
    double m_OldIntEnergy;
    Vec3D m_OldEdgeVector;    // a vector along the edge
    double m_OldEdgeSize;    // size of the edge
    int m_OldLinkType;  // 0 surface link; 1 edge link;
    std::vector<double> m_OldVFieldIntEnergy;    // interaction energy of the vector fields
    int m_Number_of_VectorField_Layers;

};


#endif
