#if !defined(AFX_Restart_H_INCLUDED_)
#define AFX_Restart_H_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This for reading and writing the restart file (a binary file)
 */
#include "SimDef.h"
#include "CreateMashBluePrint.h"

class State;
struct MESH;
class Restart
{
public:
    Restart();
	Restart(State*);
	 ~Restart();
    
public:
    
    void CopyBinaryFile(std::string file1, std::string file2); // A function to copy binary file1 into binary file2
    void CopyFile(std::string file1, std::string file2); // copy a text file into anotherone
    void WrireRestart(int step, std::string fileflag, MESH * pmesh, double r, double rb);
    MeshBluePrint ReadRestart(std::string filename, bool*);
private:
    State *m_pState;
    std::string m_CopyFileName;



};

#endif
