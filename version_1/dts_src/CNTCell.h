#if !defined(AFX_CNTCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_)
#define AFX_CNTCell_H_8P4B21B8_C13C_5648_BF23_444095086234__INCLUDED_

#include "SimDef.h"
#include "Vec3D.h"
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Unit Cell Object: A domain in space that contains a set of vertices. It is dynamic, so during simulation time, the vertices change
 It is made of an index and 3 integer to check how is staked in the space compare to other cells.
 *******************/
class vertex;
class CNTCell
{
public:
	CNTCell(int id, int i, int j, int k);
	 ~CNTCell();

// a set of functions to access the private variables
        inline const std::vector <int>   GetIndex()    const    {return m_Index;}
        inline const int GetID()                       const  	{return m_ID;}
        inline std::vector <vertex *> GetVertexList()           {return m_VertexList;} 
        inline Vec3D	GetCNTSidemin()              		    {return m_CNTRmin;}
        inline Vec3D	GetCNTSidemax()              		    {return m_CNTRmax;}
        inline std::vector <CNTCell *> GetVNeighbourCNTCell()   {return m_VNeighbourCNTCell;}

public:
    // A set of public function to change the content of the cell.
  void AddtoVertexList(vertex * z);   // adding a new vertex, when the vertex enters the domain
  void UpdateCNTSidemin(Vec3D);
  void UpdateCNTSidemax(Vec3D);
  void AddtoNeighbourCNTCell(CNTCell* z); // pointer to the nighbouring cells
  void RemoveFromVertexList(vertex * z);  // removing a  vertex, when the vertex exist the domain
   CNTCell * ReturnNewCell(int,int,int);    // accessing a specific nighbouring cell

private:
  // private member variables
  int m_ID;     // ID
  std::vector <vertex *> m_VertexList;  // List of the vertices that the domain contains
  std::vector <CNTCell *> m_VNeighbourCNTCell;  // list of the nighbouring cells
  std::vector <int>  m_Index;
  Vec3D		m_CNTRmin;
  Vec3D		m_CNTRmax;

};

#endif
