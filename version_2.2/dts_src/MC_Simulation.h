#ifndef MC_SIMULATION_H_INCLUDED
#define MC_SIMULATION_H_INCLUDED

#include "SimDef.h"
#include "AbstractSimulation.h"
/*
 * MC_Simulation Class
 * Date 2024
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 *
 * Description:
 * MC_Simulation class implements Monte Carlo simulations for the system defined by the State object.
 * It performs various simulation steps and writes output files at each step.
 *
 * Public Methods:
 * - Constructor: Initializes the MC_Simulation object with a State pointer.
 * - Destructor: Cleans up resources associated with the MC_Simulation object.
 * - Initialize: Initializes the simulation.
 * - do_Simulation: Executes the Monte Carlo simulation steps.
 *
 * Initialization Steps:
 * 1. Frame of the initial system is written.
 * 2. Mesh is voxelized for the first time.
 * 3. Mesh quality is checked.
 * 4. System curvature and energy are updated.
 *
 * Simulation Steps:
 * - Simulation progresses from the initial step to the final step.
 * - Centering of the simulation box is performed if configured.
 * - Standard Integrators are executed, including vertex position update, link flip update, and inclusion update.
 * - Supplementary Integrators, such as box update and mesh topology update, can be executed.
 * - Visualization frames, trajectory files, time series data, and checkpoint files are written based on the defined rate in the associated object.
 * - Simulation progress is printed to the console.
 *
 * Additional Notes:
 * - Energy check is performed to ensure the accuracy of calculations.
 */

class State;
class MC_Simulation : public AbstractSimulation {
public:
    MC_Simulation(State *pState);
    ~MC_Simulation();

public:
    void Initialize();
    bool do_Simulation();
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "MC_Simulation";}
    inline static std::string GetDefaultReadName() {return "MC_Simulation";}

private:
    State *m_pState;
    void PrintRate(int step, bool clean, bool clear);
    bool CheckMesh(int step);
    // Private members
};

#endif // MC_SIMULATION_H_INCLUDED
