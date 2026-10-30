#if !defined(AFX_BoxSizeCouplingToHarmonicPotentialAnisotropic_H_MM3B21B8_C13C_5648_BF23_344095086255__INCLUDED_)
#define AFX_BoxSizeCouplingToHarmonicPotentialAnisotropic_H_MM3B21B8_C13C_5648_BF23_344095086255__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractDynamicBox.h"

class State;
class BoxSizeCouplingToHarmonicPotentialAnisotropic : public AbstractDynamicBox {
public:
    
    BoxSizeCouplingToHarmonicPotentialAnisotropic(int period, double kx, double ky, double kz, double a0x, double a0y, double a0z, std::string direction, State *pState);
	~BoxSizeCouplingToHarmonicPotentialAnisotropic();

    void Initialize();
    bool ChangeBoxSize(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "HarmonicPotentialAnisotropic";}
    inline  static std::string GetDefaultReadName()  {return "HarmonicPotentialAnisotropic";}
    std::string CurrentState();

private:
    bool AnAtemptToChangeBox(double lx,double ly, double lz, double tem);
    bool VertexMoveIsFine(double lx,double ly, double lz);
    bool CheckLinkLength(double lx,double ly, double lz);
    double StretchedDistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2, double lx, double ly, double lz);
    bool CheckFaceAngleOfOneLink(links * p_edge);  // Function to check if the angle between the normal vectors of the faces sharing the given edge
    bool CheckFaces(); // Function to check if the angle between the faces of all links in the m_pRightL list
    void SetDirection(const std::string& direction);
    
private:
    Vec3D m_Direction;
    std::string m_Type;
    double m_Kx;
    double m_Ky;
    double m_Kz;
    double m_Lx0;
    double m_Ly0;
    double m_Lz0;
    Vec3D m_Box_0;
    Vec3D m_K;

    int m_Period;
    State *m_pState;
    
private:
    std::vector<vertex*>&        m_pActiveV;
    std::vector<triangle*>&      m_pActiveT;
    std::vector<links*>&   m_pRightL;
    std::vector<links*>&   m_pEdgeL;

    Vec3D *m_pBox;
    double &m_Beta;
    double &m_DBeta;
    double &m_MinLength2;
    double &m_MaxLength2;
    double &m_MinAngle;

};

/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of drx dry by calling "ChangeBoxSize" function.
=================================================================================================================
*/
#endif
