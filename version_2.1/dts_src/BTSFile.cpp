

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include <vector>
#include "BTSFile.h"
#include "State.h"
#include "Nfunction.h"
#include "CreateMashBluePrint.h"
BTSFile::BTSFile(State* pState){
    m_Periodic = 0;
    m_pState = pState;
}
BTSFile::BTSFile(State *pState, int periodic, int precision, std::string btsfilename){

    m_pState = pState;
    m_Periodic = periodic;
    m_FileName = btsfilename;
}
BTSFile::~BTSFile()
{
    m_File.close();
}
/*
bool BTSFile::OpenFile(bool clear, char readorwrite){
    
    // we could just use generic name and find the file name using state class
    std::string ext = m_FileName.substr(m_FileName.find_last_of(".") + 1);
    
    if(ext!=BTSExt)
        m_FileName = m_FileName+"."+BTSExt;
  if(clear==false && readorwrite=='w')
  m_File.open(m_FileName.c_str(),std::ios::out | std::ios::binary | std::ios::trunc );
  else if(clear==true && readorwrite=='w')
   m_File.open(m_FileName.c_str(),std::ios::out | std::ios::binary | std::ios::app );
  else if(readorwrite=='r')
    m_File.open(m_FileName.c_str(), std::ios::in |std::ios::binary);
  else{
      std::cout<<"---> error229 "<<std::endl;
      return false;
  }
    
    return true;
}*/
bool BTSFile::OpenFile(bool clear, char readorwrite) {
    // Determine file extension
    std::string ext = m_FileName.substr(m_FileName.find_last_of(".") + 1);

    // Append extension if necessary
    m_FileName = (ext != BTSExt) ? m_FileName + "." + BTSExt : m_FileName;

    // Determine file mode based on parameters
    std::ios_base::openmode mode;
    if (readorwrite == 'w') {
        mode = (clear) ? std::ios::out | std::ios::binary | std::ios::app : std::ios::out | std::ios::binary | std::ios::trunc;
    } else if (readorwrite == 'r') {
        mode = std::ios::in | std::ios::binary;
    } else {
        std::cerr << "--> error: Invalid read/write mode in bts file !" << std::endl;
        return false;
    }

    // Open the file
    m_File.open(m_FileName.c_str(), mode);

    // Check if the file opened successfully
    if (!m_File.is_open()) {
        std::cerr << "---> error: Failed to open file '" << m_FileName << "'" << std::endl;
        return false;
    }

    return true;
}
//=== Writing a BTSFile file: Writing the mesh only
void BTSFile::WriteAFrame(int &step){
    
    std::cout<<" fghgh this function is not complete yet \n";
    if(m_Periodic==0 || step%m_Periodic!=0)
        return;
        
        MESH * pmesh;
    MeshBluePrint blueprint = pmesh->Convert_Mesh_2_BluePrint(pmesh);

    (m_File).write((char *) &step, sizeof(int));
    Vec3D box = blueprint.simbox;
    (m_File).write((char *) &(box(0)), sizeof(double));
    (m_File).write((char *) &(box(1)), sizeof(double));
    (m_File).write((char *) &(box(2)), sizeof(double));

    int size = (blueprint.bvertex).size();
    (m_File).write((char *) &size, sizeof(int));
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        (m_File).write((char *) &(*it), sizeof(Vertex_Map));
    size = (blueprint.btriangle).size();
    (m_File).write((char *) &size, sizeof(int));
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        (m_File).write((char *) &(*it), sizeof(Triangle_Map));
    size = (blueprint.binclusion).size();
    (m_File).write((char *) &size, sizeof(int));
    for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        (m_File).write((char *) &(*it), sizeof(Inclusion_Map));
    
    m_File.flush();
    return;
}
//=== Read a BTSFile file
MeshBluePrint BTSFile::ReadBTSFile(bool *readok)
{
    *readok = false;
    MeshBluePrint blueprint;
#if TEST_MODE == Enabled
    std::cout<<"----> Note: BTSFile file name for reading "<<m_FileName<<std::endl;
#endif

    m_File.read((char *) &m_Step, sizeof(int));
if(m_File.is_open() && !m_File.eof() )
{
#if TEST_MODE == Enabled
    std::cout<<"----> Note: The saved BTSFile was at step "<<m_Step<<std::endl;
#endif
    double lx,ly,lz;
    m_File.read((char *) &lx, sizeof(double));
    m_File.read((char *) &ly, sizeof(double));
    m_File.read((char *) &lz, sizeof(double));
    Vec3D box(lx,ly,lz);
    blueprint.simbox= box;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector <InclusionType> binctype;  // a vector containing all inclsuion type and a default one
    int size;
    m_File.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++)
    {
        Vertex_Map vmap;
        (m_File).read((char *) &vmap, sizeof(Vertex_Map));
        bvertex.push_back(vmap);
    }
    m_File.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++)
    {
        Triangle_Map tmap;
        (m_File).read((char *) &tmap, sizeof(Triangle_Map));
        btriangle.push_back(tmap);
    }
    m_File.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++)
    {
        Inclusion_Map incmap;
        (m_File).read((char *) &incmap, sizeof(Inclusion_Map));
        binclusion.push_back(incmap);
    }

    
    blueprint.bvertex = bvertex;
    blueprint.btriangle = btriangle;
    blueprint.binclusion = binclusion;
    *readok = true;

}
    return blueprint;
}
std::string BTSFile::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}
