#if !defined(AFX_TSIFile_H_7F4B21B8_D13C_9321_QF23_885095086234__INCLUDED_)
#define AFX_TSIFile_H_7F4B21B8_D13C_9321_QF23_885095086234__INCLUDED_
#include "Def.h"
#include "Vec3D.h"
#include "Nfunction.h"
#include "MeshBluePrint.h"

class TSIFile
{
public:
    
    TSIFile();
     ~TSIFile();



public:

void WriteTSI(std::string filename, MeshBluePrint mesh);
MeshBluePrint ReadTSI(std::string filename );

private:
    std::string m_tsiPrecision;





};


#endif
