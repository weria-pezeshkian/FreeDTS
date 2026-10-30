#if !defined(AFX_GenericPositionRescaleBoxChange_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_GenericPositionRescaleBoxChange_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractDynamicBox.h"

class State;
class GenericPositionRescaleBoxChange : public AbstractDynamicBox {
public:
    
    GenericPositionRescaleBoxChange(std::istream& indata, State *pState);
	~GenericPositionRescaleBoxChange();

    void Initialize();
    bool ChangeBoxSize(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "G-PositionRescale";}
    inline  static std::string GetDefaultReadName()  {return "G-PositionRescale";}
    std::string CurrentState();

private:
    bool AnAtemptToChangeBox(double lx,double ly, double lz, double tem);
    bool VertexMoveIsFine(double lx,double ly, double lz);
    bool CheckLinkLength(double lx,double ly, double lz);
    double StretchedDistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2, double lx, double ly, double lz);
    bool CheckFaceAngleOfOneLink(links * p_edge);  // Function to check if the angle between the normal vectors of the faces sharing the given edge
    bool CheckFaces(); // Function to check if the angle between the faces of all links in the m_pRightL list
    
private:
    std::string m_Type;
    double m_SigmaP;
    int m_Period;
    State *m_pState;
    double m_Alpha;


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

private:
    bool SetInputs(std::istream& input);
    
    
    
//  Box change degree of freedom type (Semi Isotropic/)
    Vec3D (GenericPositionRescaleBoxChange::*mF_IsotropyType)();
    double (GenericPositionRescaleBoxChange::*mF_EnergyOfBoxChange)(double);

    Vec3D SemiIsotropic();
    double EnergySemiIsotropic(double lx);

};


#endif
