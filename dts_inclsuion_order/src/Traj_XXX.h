#if !defined(AFX_Traj_XXX_H_7F4B21B8_D13C_9321_QF23_885095086234__INCLUDED_)
#define AFX_Traj_XXX_H_7F4B21B8_D13C_9321_QF23_885095086234__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "Vec3D.h"
#include "inclusion.h"
#include "Nfunction.h"
#include "CreateMashBluePrint.h"

class Traj_XXX
{
public:
    
    Traj_XXX(Vec3D *pBox, std::string tsiFolder_name, std::string );
     ~Traj_XXX();





        inline int GetStep()                     			  {return m_Step;}
        inline bool GetCondition()                     			  {return m_Condition;}

public:

void WriteTSI(int step ,  std::string filename, std::vector< vertex* > pver, std::vector< triangle* > ptriangle,   std::vector< inclusion* > pinc);
MeshBluePrint ReadTSI1(std::string filename);
MeshBluePrint ReadTSI2(std::string filename, std::vector<InclusionType> bINCtype);
MeshBluePrint ReadTSI(std::string filename , std::vector<InclusionType> bINCtype, bool*);


private:


    bool m_Condition;
    Vec3D *m_pBox;
    std::string m_filename;
    int m_Step;
    std::string m_tsiFolder_name;
    std::string m_tsiPrecision;




};


#endif
