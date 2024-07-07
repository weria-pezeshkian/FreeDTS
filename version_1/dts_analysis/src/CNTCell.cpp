

#include <stdio.h>
#include "CNTCell.h"
CNTCell::CNTCell(int id, int i, int j, int k)
{
m_Index.push_back(i);
m_Index.push_back(j);
m_Index.push_back(k);
m_ID=id;

}
CNTCell::~CNTCell()
{
    
}
void CNTCell::AddtoNeighbourCNTCell(CNTCell* z)
{
m_VNeighbourCNTCell.push_back(z);
}
void CNTCell::AddtoVertexList(vertex * z)
{

m_VertexList.push_back(z);
}
void CNTCell::UpdateCNTSidemin(Vec3D z)
{

m_CNTRmin=z;
}
void CNTCell::UpdateCNTSidemax(Vec3D z)
{

m_CNTRmax=z;
}

void CNTCell::RemoveFromVertexList(vertex * z)
{
m_VertexList.erase(std::remove(m_VertexList.begin(), m_VertexList.end(), z), m_VertexList.end());
}
CNTCell * CNTCell::ReturnNewCell(int i,int j,int k)
{
    int n=0;
    n=(i+1)+3*(j+1)+9*(k+1);
    if(n>13)
    {
        n=n-1;
    }
    CNTCell * P=m_VNeighbourCNTCell.at(n);
    return P;
}



