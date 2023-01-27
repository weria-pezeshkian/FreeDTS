#if !defined(AFX_GenerateCNTCells_H_8P4B21B8_C13C_5648_BF23_444095086239__INCLUDED_)
#define AFX_GenerateCNTCells_H_8P4B21B8_C13C_5648_BF23_444095086239__INCLUDED_


#include "SimDef.h"
#include "CNTCell.h"
#include "vertex.h"
#include "Vec3D.h"
class GenerateCNTCells
{
public:
    

	GenerateCNTCells(std::vector< vertex* > ver,Vec3D cntsize,Vec3D *pBox);
	~GenerateCNTCells();




        inline std::vector <CNTCell *> GetAllCNTCells()           {return m_pAllCNTCells;}
    inline std::vector <double> GetCNTCellSize()        {return m_CNTCellSize;}
    inline std::vector <int> GetCNTCellNo()        {return m_CNTCellNo;}




public:
    

//  void AddtoVertexList(vertex * z);

    void Generate();

private:
std::vector <CNTCell *> m_pAllCNTCells;
std::vector <CNTCell > m_AllCNTCells;
int IDFromIndex(int,int,int);
int IndexFromID(int,int *,int *,int *);
std::vector <double> m_CNTSize;
private:
    std::vector< vertex* > m_Allver;

int m_Nx;
int m_Ny;
int m_Nz;
    Vec3D  *m_pBox;
    std::vector <double> m_CNTCellSize;
    std::vector <int> m_CNTCellNo;




    





};


#endif
