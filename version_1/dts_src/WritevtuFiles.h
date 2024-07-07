#if !defined(AFX_WritevtuFiles_H_7F4B21B8_D13C_9321_BF23_124095086234__INCLUDED_)
#define AFX_WritevtuFiles_H_7F4B21B8_D13C_9321_BF23_124095086234__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "Vec3D.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An object for writing vtu files (to visualize the meshes using paraview).
 */
class WritevtuFiles
{
public:
    
	  WritevtuFiles(State* pState);
	 ~WritevtuFiles();
public:

void Writevtu(std::vector<vertex* > ver, std::vector<triangle* > triangle,  std::vector<links* > , std::string Filename);
    void Writevtu_Plus_NumberList(std::vector<std::string> num_name, std::vector<std::vector<double> > num, std::vector<vertex* > ver, std::vector<triangle* > triangle,  std::vector<links* > , std::string Filename);  // using this, we can write a vtu file that has energy list or anything 
void WritevtuNochange(std::vector<vertex* > ver, std::vector<triangle* > triangle,  std::vector<links* > , std::string Filename);
void Writefullvtu(std::vector<vertex* > ver, std::vector<triangle* > triangle,  std::vector<links* > , std::string Filename);
void WriteInclusion(std::string id, std::vector<vertex* > ver, std::ofstream *Output);
private:

Vec3D *m_pBox;
    State* m_pState;

};
#endif
