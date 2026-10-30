#ifndef RESTART_H_INCLUDED
#define RESTART_H_INCLUDED

#include "SimDef.h"
#include "CreateMashBluePrint.h"
/*
 * @brief Header file for the Restart class, which manages reading and writing of a restart file for simulation states.
 *
 * This header file contains the declaration of the Restart class and its associated data structures and methods.
 * The Restart class is responsible for saving the current state of a simulation to a binary restart file and
 * reading from such a file to resume the simulation.
 *
 * @details The Restart class provides functionality to:
 * - Update the restart state with the current simulation step and parameters.
 * - Read the simulation state from a restart file.
 * - private function: Convert between string representations of vector fields and their structured data forms.
 * - private function: Set and manage the names of temporary and final restart files.
 *
 * The header file also defines the SingleVF structure used to represent vector fields in the simulation.
 *
 * @author
 * Weria Pezeshkian (weria.pezeshkian@gmail.com)
 *
 * @copyright
 * Copyright (c) Weria Pezeshkian
 */
class State;
struct MESH;

struct SingleVF {
    // Data structure for vector field map local to this class. This is not needed anywhere else
    // it is different from vector field map in the blue print
    // The SingleVF structure is used to store information about a single vector field, including
    // its inclusion type and its x and y components.
    int inc_type; // The inclusion type of the vector field.
    double x; // The x-component of the vector field.
    double y; // The y-component of the vector field.
};
class Restart {
public:
    Restart();
    Restart(State* pstate);
    ~Restart();
    
    
    bool UpdateRestartState(int step, double r_vertex, double r_box);  // Writing the restart file
    /*
     * @brief Updates the restart state by writing the current simulation step and parameters to a file.
     *
     * @param step The current simulation step.
     * @param r_vertex A double value representing a parameter related to vertices.
     * @param r_box A double value representing a parameter related to the simulation box.
     * @return True if the update is successful, false otherwise.
     */
    MeshBluePrint ReadFromRestart(const std::string& filename, int& step, bool& readok, double& r_vertex, double& r_box); // Reading a restart file
    void SetRestartFileName();
    void UpdatePeriod(int period);
    std::string CurrentState();

private:
    void WriteRestart(std::string &filename, int step, MESH* pmesh, double r, double rb);
    MeshBluePrint ReadRestart(std::string filename, int& step, bool& readok, double& r_vertex, double& r_box); // Reading a restart file

    std::vector<SingleVF> convert2VF(int no_vf, std::string data);   // brief: Converts a string representation of vector fields into a vector of SingleVF objects.
    std::string V_VF2String(std::vector<SingleVF>  v_VF); // brief: Converts a vector of SingleVF objects into a string representation.

    State* m_pState;
    std::string m_TEMFileName;
    std::string m_RestartFileName;
    int m_Period;
};

#endif // RESTART_H_INCLUDED
