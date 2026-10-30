#if !defined(AFX_BTSFile_H_INCLUDED_)
#define AFX_BTSFile_H_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This for reading and writing the BTSFile file (a binary file for trajectory)
 */
#include "SimDef.h"
#include "CreateMashBluePrint.h"
#include "AbstractBinaryTrajectory.h"

class State;
class BTSFile : public  AbstractBinaryTrajectory {
public:
    BTSFile(State *pstate);
    BTSFile(State *pstate, int periodic, int precision, std::string btsfilename);
    
	 ~BTSFile();
    
public:
    bool OpenFile(bool clear, char rw);
    void WriteAFrame(int &step);
    MeshBluePrint ReadBTSFile(bool*);
    
    inline  std::string GetDerivedDefaultReadName() {return "OutPutTRJ_BTS";}
    inline static std::string GetDefaultReadName() {return "OutPutTRJ_BTS";}
    std::string CurrentState();

private:
    int m_Step;
    std::string m_FileName;
    std::fstream m_File;
    State *m_pState;
    int m_Periodic;

};

#endif
