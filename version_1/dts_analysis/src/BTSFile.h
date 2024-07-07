#if !defined(AFX_BTSFile_H_INCLUDED_)
#define AFX_BTSFile_H_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This for reading and writing the BTSFile file (a binary file for trajectory)
 */
#include "SimDef.h"
#include "CreateMashBluePrint.h"

class MESH;
class BTSFile
{
public:
    BTSFile();
    BTSFile(std::string filename, bool clear, std::string writeorread);

	 ~BTSFile();
    
public:    
    void WrireBTSFile(int step, MESH * pmesh);
    MeshBluePrint ReadBTSFile(bool*);
private:
    int m_Step;
    std::string m_FileName;
    std::fstream m_File;


};

#endif
