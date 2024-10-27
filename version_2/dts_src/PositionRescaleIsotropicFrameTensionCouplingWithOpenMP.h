#if !defined(AFX_PositionRescaleIsotropicFrameTensionCouplingWithOpenMP_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_PositionRescaleIsotropicFrameTensionCouplingWithOpenMP_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractDynamicBox.h"

class State;
class PositionRescaleIsotropicFrameTensionCouplingWithOpenMP : public AbstractDynamicBox {
public:
    
    PositionRescaleIsotropicFrameTensionCouplingWithOpenMP(int period, double tau, std::string direction, State *pState);
	~PositionRescaleIsotropicFrameTensionCouplingWithOpenMP();

    void Initialize();
    bool ChangeBoxSize(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "IsotropicFrameTensionOpenMP";}
    inline  static std::string GetDefaultReadName()  {return "IsotropicFrameTensionOpenMP";}
    std::string CurrentState();

private:
    bool AnAtemptToChangeBox(double lx,double ly, double lz, double tem);
    bool VertexMoveIsFine(double lx,double ly, double lz);
    bool CheckLinkLength(double lx,double ly, double lz);
    double StretchedDistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2, double lx, double ly, double lz);
    bool CheckFaceAngleOfOneLink(links * p_edge);  // Function to check if the angle between the normal vectors of the faces sharing the given edge
    bool CheckFaces(); // Function to check if the angle between the faces of all links in the m_pRightL list
    void SetDirection(std::string direction);
    
private:
    Vec3D m_Direction;
    std::string m_Type;
    double m_SigmaP;
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


#endif
