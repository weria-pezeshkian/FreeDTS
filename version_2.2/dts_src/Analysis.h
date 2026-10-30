#ifndef Analysis_H_INCLUDED
#define Analysis_H_INCLUDED

#include "SimDef.h"
#include "AbstractSimulation.h"
/*
 * Analysis Class
 * Date 2024
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 *
 * Description:
 * Analysis class implements Monte Carlo simulations for the system defined by the State object.
 * It performs various simulation steps and writes output files at each step.
 *
 * Public Methods:
 * - Constructor: Initializes the Analysis object with a State pointer.
 * - Destructor: Cleans up resources associated with the Analysis object.
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
class Analysis : public AbstractSimulation {
public:
    Analysis(State *pState, std::string path);
    ~Analysis();

public:
    void Initialize();
    bool do_Simulation();
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName() {return "Analysis";}
    inline static std::string GetDefaultReadName() {return "Analysis";}

private:
    State *m_pState;
    bool CheckMesh(int step);
    std::string m_Path;
    // Private members
    
private:
    bool CenterMesh();
};

#endif // Analysis_H_INCLUDED
