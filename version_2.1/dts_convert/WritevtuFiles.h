#if !defined(AFX_WritevtuFiles_H_7F4B21B8_D13C_9321_BF23_124095086234__INCLUDED_)
#define AFX_WritevtuFiles_H_7F4B21B8_D13C_9321_BF23_124095086234__INCLUDED_

#include "Def.h"
#include "Vec3D.h"
#include "MeshBluePrint.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An object for writing vtu files (to visualize the meshes using paraview).
 */
class WritevtuFiles
{
public:
    
	  WritevtuFiles();
	 ~WritevtuFiles();
public:

    void Write(MeshBluePrint blue, std::string filename);
    bool CrossPBC(Vec3D Box, Vertex_Map v1, Vertex_Map v2, Vertex_Map v3);
private:

};
#endif
