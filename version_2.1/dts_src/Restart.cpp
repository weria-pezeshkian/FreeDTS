

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Restart.h"
#include "State.h"
#include "CreateMashBluePrint.h"

Restart::Restart(){
    m_Period = 1000;
    m_RestartFileName = "";
}
Restart::Restart(State *pState)
{
    m_pState = pState;
    m_TEMFileName = "-1.res";
    m_Period = 1000;
    m_RestartFileName = ".res";
}
Restart::~Restart()
{
    
}
void Restart::UpdatePeriod(int period){
    
    m_Period = period;
    return;
}
void Restart::SetRestartFileName(){
    
    m_RestartFileName = m_pState->GetRunTag() +"."+ RestartExt;
    m_TEMFileName = m_pState->GetRunTag() +"-1."+ RestartExt;
    return;
}
bool Restart::UpdateRestartState(int step,  double r_vertex, double r_box){
        
    if(m_Period!=0 && step%m_Period == 0 ){
        WriteRestart(m_RestartFileName, step, m_pState->GetMesh(),r_vertex, r_box);
    }
    return true;
}
//=== Writing a restart file: this file is just the active [State] Object with an updated initial step
void Restart::WriteRestart(std::string &filename, int step, MESH * pmesh, double r, double rb){
    
    /**
     * @brief Writes a restart file containing the current state of the simulation (it need to be coupled to right input file).
     *
     * This function writes the state of the simulation into a binary file, including the current
     * step, mesh data, simulation parameters, and vector fields. This restart file allows the
     * simulation to be resumed from the saved state.
     *
     * @param filename The name of the restart file to be written.
     * @param step The current step of the simulation.
     * @param pmesh Pointer to the MESH object representing the mesh data.
     * @param r A double value representing a simulation parameter.
     * @param rb A double value representing another simulation parameter.
     *
     * Example usage:
     * @code
     * Restart restart;
     * std::string filename = "restart.res";
     * int step = 100;
     * MESH mesh;
     * double r = 1.0;
     * double rb = 0.5;
     * restart.WriteRestart(filename, step, &mesh, r, rb);
     * @endcode
     */

    MeshBluePrint blueprint = pmesh->Convert_Mesh_2_BluePrint(pmesh);

    //=== we write the restart into tem restart file
    std::fstream Rfile;
    Rfile.open(m_TEMFileName.c_str(),std::ios::out | std::ios::binary);
    (Rfile).write((char *) &step, sizeof(int));    
    (Rfile).write((char *) &r, sizeof(double));
    (Rfile).write((char *) &rb, sizeof(double));
    Vec3D box = blueprint.simbox;
    (Rfile).write((char *) &(box(0)), sizeof(double));
    (Rfile).write((char *) &(box(1)), sizeof(double));
    (Rfile).write((char *) &(box(2)), sizeof(double));

    int size = (blueprint.bvertex).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it){
        // write vertices
        (Rfile).write((char *) &(*it), sizeof(Vertex_Map));
    }
    size = (blueprint.btriangle).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        (Rfile).write((char *) &(*it), sizeof(Triangle_Map));
    size = (blueprint.binclusion).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        (Rfile).write((char *) &(*it), sizeof(Inclusion_Map));
    
    //-- write vector fields
    int no_vec = blueprint.number_vector_field;
    (Rfile).write((char *) &no_vec, sizeof(int));
    
    for (std::vector<VectorField_Map>::iterator it = (blueprint.bvectorfields).begin() ; it != (blueprint.bvectorfields).end(); ++it){
        std::vector<SingleVF> V_vf = convert2VF(no_vec, it->data_line);
        for (int i = 0; i<no_vec; i++){
            SingleVF tm_vf = V_vf[i];
            (Rfile).write((char *) &tm_vf, sizeof(SingleVF));
        }
    }
    // Close the temporary restart file and flush the buffer
    Rfile.close();
    Rfile.flush();
    // Copy the temporary restart file to the final restart file
    Nfunction::CopyBinaryFile(m_TEMFileName,m_RestartFileName, FreeDTS_BUFFERSIZE);
    // Remove the temporary restart file
    remove(m_TEMFileName.c_str());

    return;
}
MeshBluePrint Restart::ReadFromRestart(const std::string& inputFilename, int& step, bool& readOk, double& rVertex, double& rBox) {
    readOk = false;
    std::string restartFilename = m_TEMFileName; // Default to temporary file name

    // Check if the temporary restart file exists
    std::ifstream tempFile(m_TEMFileName);
    if (!tempFile.good()) {
        // Use the provided input filename if the temporary file doesn't exist
        restartFilename = inputFilename;

        // Append the restart file extension if needed
        std::string ext = inputFilename.substr(inputFilename.find_last_of(".") + 1);
        if (ext != RestartExt) {
            restartFilename = restartFilename + "." + RestartExt;
        }
    }

    // Read from the determined restart file
    return ReadRestart(restartFilename, step, readOk, rVertex, rBox);
}
//=== Read a restart file and load to the  active [State] Object
MeshBluePrint Restart::ReadRestart(std::string filename , int &step, bool &readok, double &r_vertex, double &r_box) {

    MeshBluePrint blueprint;

    // We should make a check and see if the restart file is writtn correctly or not
    //=== just read the restart file and copy it into the current active [State] object
    std::fstream Rfile;
    Rfile.open(filename.c_str(), std::ios::in |std::ios::binary);
if(Rfile.is_open())
{
    Rfile.read((char *) &step, sizeof(int));

    double lx,ly,lz;
    Rfile.read((char *) &r_vertex, sizeof(double));
    Rfile.read((char *) &r_box, sizeof(double));
    Rfile.read((char *) &lx, sizeof(double));
    Rfile.read((char *) &ly, sizeof(double));
    Rfile.read((char *) &lz, sizeof(double));
    Vec3D box(lx,ly,lz);
    blueprint.simbox= box;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    std::vector<VectorField_Map> bvectorfields; // a vector of all inclusions (only the blueprint not the object) in the mesh

    int no_ver;
    Rfile.read((char *) &no_ver, sizeof(int));
    for (int i=0;i<no_ver;i++){
        Vertex_Map vmap;
        (Rfile).read((char *) &vmap, sizeof(Vertex_Map));
        vmap.include = true;   // 
        bvertex.push_back(vmap);
    }
    int size;
    Rfile.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++){
        Triangle_Map tmap;
        (Rfile).read((char *) &tmap, sizeof(Triangle_Map));
        btriangle.push_back(tmap);
    }
    Rfile.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++){
        Inclusion_Map incmap;
        (Rfile).read((char *) &incmap, sizeof(Inclusion_Map));
        binclusion.push_back(incmap);
    }
    //-- read vector fields
    int no_vec;
    Rfile.read((char *) &no_vec, sizeof(int));
    
    for (int nv = 0; nv < no_ver; nv++){
        std::vector<SingleVF> V_vf;
        for (int i = 0; i<no_vec; i++){
            SingleVF tm_vf;
            (Rfile).read((char *) &tm_vf, sizeof(SingleVF));
            V_vf.push_back(tm_vf);
        }
        VectorField_Map t_vfm;
        t_vfm.data_line = V_VF2String(V_vf);
        bvectorfields.push_back(t_vfm);
    }

    blueprint.bvertex = bvertex;
    blueprint.btriangle = btriangle;
    blueprint.binclusion = binclusion;
    blueprint.number_vector_field = no_vec;
    blueprint.bvectorfields = bvectorfields;

    Rfile.close();

}
    readok = true;
    return blueprint;
}
std::vector<SingleVF> Restart::convert2VF(int no_vfm, std::string data){   //
    /**
     * @brief Converts a string representation of vector fields into a vector of SingleVF objects.
     *
     * This function takes a string containing vector field data and converts it into a vector
     * of SingleVF objects. The string data is expected to be formatted such that each vector
     * field is represented by three consecutive substrings: the inclusion type, the x-coordinate,
     * and the y-coordinate. The number of vector fields (no_vfm) is also provided to ensure the
     * correct parsing of the string.
     *
     * @param no_vfm The number of vector fields to be converted.
     * @param data A string containing the vector field data.
     * @return std::vector<SingleVF> A vector of SingleVF objects representing the vector fields.
     *         If the input data does not match the expected format, an empty vector is returned.
     *
     * @note The input data string is split into substrings using the Nfunction::Split method.
     *       Each vector field is then constructed from three consecutive substrings:
     *       the inclusion type, the x-coordinate, and the y-coordinate.
     *       If the number of substrings does not match the expected 3 * no_vfm, an error message
     *       is printed and an empty vector is returned.
     *
     * Example usage:
     * @code
     * Restart restart;
     * std::string data = "1 0.5 0.6 2 0.7 0.8";
     * int no_vfm = 2;
     * std::vector<SingleVF> vectorFields = restart.convert2VF(no_vfm, data);
     * @endcode
     */
    std::vector<SingleVF> VF;
    std::vector<std::string> data_str = Nfunction::Split(data);

    // Check if the data provided matches the expected format.
    if (data_str.size() != 3 * no_vfm) {
        std::cout << "---> error (restart), info on vector field is not enough \n";
    }

    // Create.
    for (int i = 0; i < no_vfm; ++i) {

        SingleVF t_vf;
        t_vf.inc_type = Nfunction::String_to_Int(data_str[i * 3]);
        t_vf.x = Nfunction::String_to_Double(data_str[i * 3 + 1]);
        t_vf.y = Nfunction::String_to_Double(data_str[i * 3 + 2]);
        VF.push_back(t_vf);
    }
    
    
    return VF;
}
std::string Restart::V_VF2String(std::vector<SingleVF>  v_VF){
    /**
     * @brief Converts a vector of SingleVF objects into a string representation.
     *
     * This function takes a vector of SingleVF objects and converts it into a string representation.
     * Each SingleVF object contains three components: the inclusion type, the x-coordinate, and the
     * y-coordinate. The resulting string contains these components separated by spaces.
     *
     * @param v_VF A vector of SingleVF objects to be converted.
     * @return std::string A string representation of the vector of SingleVF objects.
     *
     * @note The output string contains the inclusion type, x-coordinate, and y-coordinate of each
     *       SingleVF object, separated by spaces. Each set of three components is also separated by a space.
     *
     * Example usage:
     * @code
     * Restart restart;
     * std::vector<SingleVF> vectorFields = { {1, 0.5, 0.6}, {2, 0.7, 0.8} };
     * std::string data = restart.V_VF2String(vectorFields);
     * @endcode
     */
    std::string out_str;
    
    for (std::vector<SingleVF>::iterator it = v_VF.begin() ; it != v_VF.end(); ++it){
        
        out_str += Nfunction::D2S(it->inc_type) +" "+Nfunction::D2S(it->x) +" "+Nfunction::D2S(it->y) +" ";
    }
    return out_str;
}
std::string Restart::CurrentState(){
        std::string state = "Restart_Period  =  "+Nfunction::D2S(m_Period);
        return state;
}
