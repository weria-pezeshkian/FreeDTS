#ifndef TRAJ_TSI_H_INCLUDED
#define TRAJ_TSI_H_INCLUDED

#include "SimDef.h"
#include "AbstractNonbinaryTrajectory.h"
/**
 * this class writes the state of the system at a specific simulation step to a TSI file.
 * This function writes the current state of the system, including vertex positions,
 * triangle connectivity, and inclusion information, to a TSI (Triangulated Surface + 
 * Interaction) file format. The function takes the simulation step number and the
 * desired filename as input parameters. It checks if the provided step number
 * aligns with the specified write period before proceeding to write the frame.
 *
 * If the provided filename does not have the TSI file extension, ".tsi" is appended
 * to the filename to ensure proper file naming convention.
 */
// Forward declaration
class State;

class Traj_tsi : public AbstractNonbinaryTrajectory {
public:
    // Constructors and Destructor
    Traj_tsi(State *pstate);
    Traj_tsi(State *pstate, int period, std::string tsiFolder_name);
    ~Traj_tsi();

    // Public member functions
    void WriteAFrame(std::string filename);
    void WriteAFrame(int step);
    bool OpenFolder();
    
    inline  std::string GetDerivedDefaultReadName()  {return "TSI";}
    static inline  std::string GetDefaultReadName()  {return "TSI";}
    std::string CurrentState();

private:
    // Private member variables
    std::string m_Folder_name;
    std::string m_Precision;
    int m_Period;
    State* m_pState;

private:  // Consider removing this section if not needed
    // Private member functions
   // MeshBluePrint ReadAFrame(std::string filename, bool& readok);
   // MeshBluePrint ReadTSI2(std::string filename, std::vector<InclusionType> bINCtype);
};

#endif // TRAJ_TSI_H_INCLUDED
