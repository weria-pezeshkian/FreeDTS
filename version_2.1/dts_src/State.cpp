
#include <fstream>
#include <time.h>
#include "MESH.h"
#include "CreateMashBluePrint.h"
#include "State.h"
State::State(){
    
}
State::State(std::vector<std::string> argument) :

      m_NumberOfErrors(0),
      m_NumberOfWarnings(0),
      //============ Initialization of all files
// make voxels
      m_pVoxelization(new Voxelization<vertex>()),
      m_pSimulation(new MC_Simulation(this)),
      m_pNonequilibriumCommands(new NonequilibriumCommands(this)),
      m_pTimeSeriesDataOutput(new TimeSeriesDataOutput(this)),  // Initialize TimeSeriesDataOutput
      m_pTimeSeriesLogInformation(new TimeSeriesLogInformation(this)),  // Initialize TimeSeriesLogInformation
      m_pRestart(new Restart(this)),  // Initialize Restart
      m_pNonbinaryTrajectory(new Traj_tsi(this)),  // Initialize Traj_tsi
      m_pBinaryTrajectory(new NoFile),  // Initialize BinaryTrajectory

    //---- nonbonded interactions
    m_NonbondedInteractionBetweenVertices(new NoNonbondedInteractionBetweenVertices),
      //---- Initialize constraint components
      m_pVAHCalculator(new VAHGlobalMeshProperties),
      m_pVolumeCoupling(new NoCoupling(m_pVAHCalculator)),
      m_pCoupleGlobalCurvature(new NoGlobalCurvature(m_pVAHCalculator)),
      m_pTotalAreaCoupling(new NoTotalAreaCoupling(m_pVAHCalculator)),
      m_pForceonVerticesfromInclusions(new NoForce),  // Initialize ForceonVerticesfromInclusions
      m_pForceonVertices(new NoForceonVertices),
      m_pForceonVerticesfromVectorFields(new NoVFForce),  // Initialize ForceonVerticesfromvectorfields
      m_pExternalFieldOnVectorFields(new NoExternalFieldOnVectorFields),  // Initialize ExternalFieldOnVectorFields
      m_pExternalFieldOnInclusions(new NoExternalFieldOnInclusions),  // Initialize ExternalFieldOnVectorFields
      m_pVertexAdhesionToSubstrate(new NoVertexAdhesionCoupling),  // Initialize ExternalFieldOnVectorFields
      m_pApplyConstraintBetweenGroups(new NoConstraint),  // Initialize ApplyConstraintBetweenGroups
      m_pBoundary(new PBCBoundary),  // Initialize Boundary
      m_pBondedPotentialBetweenVertices(new EmptyBonds),

//---- Initialize local stretching
      m_pAbstractLocalStretching(new NoLocalStretching),
      //---- Initialize supplementary integrators
      m_pDynamicBox(new NoBoxChange),  // Initialize DynamicBox
      m_pDynamicTopology(new ConstantTopology),  // Initialize DynamicTopology
      m_pOpenEdgeEvolution(new NoEvolution),  // Initialize OpenEdgeEvolution
      m_pInclusionConversion(new NoInclusionConversion),  // Initialize InclusionConversion

      //--- Initialize accessory objects
      m_pCurvatureCalculations(new CurvatureByShapeOperatorType1(this)),  // Initialize CurvatureCalculations
      m_RandomNumberGenerator(new RNG(1234)),  // Initialize RandomNumberGenerator

        // Initialize Simulation
     //--- Local state member variables
      m_Argument(argument),  // Set Argument
      m_TopologyFile("topology.top"),  // Set TopologyFile
      m_InputFileName("Input.dts"),  // Set InputFileName
      m_IndexFileName("Index.inx"),  // Set IndexFileName
      m_GeneralOutputFilename("dts"),  // Set GeneralOutputFilename
      m_Total_no_Threads(1),  // Set Total_no_Threads
      m_Targeted_State(true),  // Set Targeted_State
      m_RestartFileName("")  // Set RestartFileName
{ // start of the "class State" constructor
    
#ifdef _OPENMP
    m_Total_no_Threads = omp_get_num_procs(),
    omp_set_num_threads(m_Total_no_Threads),
#endif
    //---- Initialize integrators
    m_pMesh = &m_Mesh;
    m_pEnergyCalculator         = new Energy(this);  // Initialize EnergyCalculator
    m_pVisualizationFile        = new WritevtuFiles(this);  // Initialize WritevtuFiles
    m_pVertexPositionIntegrator = new EvolveVerticesByMetropolisAlgorithm(this);
    m_pAlexanderMove            = new AlexanderMoveByMetropolisAlgorithm(this);         // Initialize AlexanderMove
    m_pInclusionPoseIntegrator  = new InclusionPoseUpdateByMetropolisAlgorithm(this);  // Initialize InclusionPoseIntegrator
    m_pVectorFieldsRotationIntegrator  = new VectorFieldsRotationByMetropolisAlgorithm (this);  // Initialize InclusionPoseIntegrator

    m_Parallel_Replica.State    = false;

//-- Find input files name (input.tsi, top.top||top.tsi, index.inx, restart.res)
    if (!ExploreArguments(argument)) {
            exit(0);
     }
     if (!ReadInputFile(m_InputFileName)) {
            exit(0);
     }
     else{
         
         ExploreArguments(argument);      // Update last state
     }
}
State::~State()
{
    delete m_pTimeSeriesDataOutput;
    delete m_pTimeSeriesLogInformation;
    delete m_pRestart;
    delete m_pNonbinaryTrajectory;
    delete m_pVisualizationFile;
    delete m_pBinaryTrajectory;
    delete m_pVAHCalculator;
    delete m_pVolumeCoupling;
    delete m_pCoupleGlobalCurvature;
    delete m_pTotalAreaCoupling;
    delete m_pForceonVerticesfromInclusions;
    delete m_pForceonVerticesfromVectorFields;
    delete m_pForceonVertices;
    delete m_pExternalFieldOnVectorFields;
    delete m_pExternalFieldOnInclusions;
    delete m_pVertexAdhesionToSubstrate;
    delete m_pApplyConstraintBetweenGroups;
    delete m_pBoundary;
    delete m_pVertexPositionIntegrator;
    delete m_pAlexanderMove;
    delete m_pInclusionPoseIntegrator;
    delete m_pVectorFieldsRotationIntegrator;
    delete m_pDynamicBox;
    delete m_pDynamicTopology;
    delete m_pOpenEdgeEvolution;
    delete m_pInclusionConversion;
    delete m_RandomNumberGenerator;
    delete m_pEnergyCalculator;
    delete m_pSimulation;
    delete m_pVoxelization;
    delete m_pNonequilibriumCommands;
    delete m_pAbstractLocalStretching;
}
bool State::ExploreArguments(std::vector<std::string> &argument){
    /*
     * ExploreArguments Function
     * Description:
     * This function parses the command-line arguments passed to the program and updates the State object accordingly.
     * The function iterates through the vector of arguments, processing each flag and its corresponding value.
     *
     * Arguments:
     * - argument: A vector of strings containing command-line arguments and their values.
     *
     * Flags and Corresponding Actions:
     * - "-h" or "--help": Displays a help message and returns false.
     * - "-in": Sets the input file name in the State object.
     * - "-top": Sets the topology file name in the State object.
     * - "-defout": Sets the general output filename in the State object.
     * - "-b": Updates the initial step value in the simulation object.
     * - "-e": Updates the final step value in the simulation object.
     * - "-seed": Updates the seed value in the simulation object.
     * - "-restart": Sets the restart filename in the State object.
     * - "-nt": Updates the total number of threads in the State object.
     * - "-ndx": Sets the index filename in the State object.
     *
     * Return Value:
     * - True if all arguments are successfully processed, false otherwise.
     *
     * Additional Notes:
     * - The function checks for unknown command-line arguments and displays an error message.
     * - It provides a usage message if an unknown flag is encountered.
     * - String_to_Int function is used to convert string arguments to integer values.
     */
    const std::string HELP_FLAG        = "-h";
    const std::string LONG_HELP_FLAG   = "--help";
    const std::string INPUT_FLAG       = "-in";
    const std::string TOPOLOGY_FLAG    = "-top";
    const std::string INI_STEP_FLAG    = "-b";
    const std::string FINAL_STEP_FLAG  = "-e";
    const std::string SEED_FLAG        = "-seed";
    const std::string RESTART_FLAG     = "-restart";
    const std::string INDEX_FLAG       = "-ndx";
    const std::string DEFOUT_FLAG      = "-defout";
    const std::string THREAD_FLAG      = "-nt";
    const std::string ANALYSIS_FLAG    = "-analysis";

    for (size_t i=1;i<argument.size();i=i+2)
    {
        const std::string& flag = argument[i];
        
        if (flag == HELP_FLAG || flag == LONG_HELP_FLAG) {
            Nfunction::HelpMessage(); // write a help message
            return false;
        }
        else if(flag == INPUT_FLAG) {
            
            m_InputFileName = argument[i+1];
        }
        else if(flag == TOPOLOGY_FLAG) {
            
            m_TopologyFile = argument[i+1];
        }
        else if(flag == DEFOUT_FLAG)
        {
            m_GeneralOutputFilename = argument[i+1];
        }
        else if(flag == INI_STEP_FLAG){
            
            m_pSimulation->UpdateInitialStep(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == FINAL_STEP_FLAG){
            
            m_pSimulation->UpdateFinalStep(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == SEED_FLAG){
            
            m_RandomNumberGenerator = new RNG(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == ANALYSIS_FLAG){
            m_pSimulation = new Analysis(this,argument[i+1] );
        }
        else if(flag == RESTART_FLAG){
            
            m_RestartFileName = argument[i+1];
        }
        else if(flag == THREAD_FLAG){
            
            m_Total_no_Threads = Nfunction::String_to_Int(argument[i+1]);
            #ifdef _OPENMP
            if(omp_get_num_procs()<m_Total_no_Threads){
                std::cout<<"---> warning, the number of requested threads is larger then the total physical cores \n";
            }
            omp_set_num_threads(m_Total_no_Threads); // Set the number of threads
            #endif
            //std::cout<<omp_get_num_procs()<<"\n";

        }
        else if(flag == INDEX_FLAG){   // get the index file and also check if the file exist
            m_IndexFileName = argument[i+1];
        }
        else
        {
            std::cout << "--> error: unknown commandline argument! "<<flag<<std::endl;
            std::cout << "    For more information and tips run "<< EXE_NAME <<" -h"<<std::endl;
            return false;
        }
    }

    return true;
}
bool State::ReadInputFile(std::string file)
{
    /*
     * @brief Reads and parses the input file for the State class, configuring various simulation parameters.
     *
     * This function reads an input file specified by the `file` parameter, which is expected to contain a variety
     * of simulation parameters and settings. The file format is assumed to follow a specific structure, where each
     * line starts with a keyword that determines the type of data or configuration being set. The function handles
     * various keywords to initialize different aspects of the simulation, such as simulation types, boundary conditions,
     * integrators, and other parameters.
     *
     * The function performs the following main steps:
     *
     * 1. **File Extension Check**:
     *    - Checks if the provided file name has the correct extension. If not, appends the expected extension.
     *
     * 2. **File Existence Check**:
     *    - Verifies if the file exists. If not, outputs an error message and returns `false`.
     *
     * 3. **File Opening**:
     *    - Opens the file for reading. If the file cannot be opened, outputs an error message and returns `false`.
     *
     * 4. **Reading and Parsing the File**:
     *    - Reads the file line by line, using the first word of each line (keyword) to determine the type of data
     *      and performs the appropriate configuration.
     *    - For each recognized keyword, extracts the relevant parameters and initializes or updates the corresponding
     *      class members or settings.
     *    - Handles various types of configuration blocks, such as:
     *      - Simulation settings (e.g., Integrator type, steps, temperature).
     *      - Visualization settings.
     *      - Integrator settings for vertex positions.
     *      - Force fields and constraints.
     *      - Boundary conditions.
     *      - Topology and dynamic box settings.
     *      - Various other simulation parameters (e.g., curvature, volume constraints, area coupling).
     *
     * 5. **Inclusion Handling**:
     *    - Calls the `ReadInclusionType` function to handle inclusion types specified in the file.
     *
     * 6. **File Closing**:
     *    - Closes the input file.
     *
     * 7. **Return Value**:
     *    - Returns `true` if the file is successfully read and parsed without errors.
     *    - Returns `false` if any errors occur during the reading or parsing process.
     *
     * @param file The name of the input file to be read.
     * @return `true` if the file is successfully read and parsed, `false` otherwise.
     */
    std::string ext = file.substr(file.find_last_of(".") + 1);
    std::string filename = (ext != InExt) ? file + "." + InExt : file;

    if (!Nfunction::FileExist(filename)) {
        std::cerr << "----> Error: the input file with the name " << filename << " does not exist" << std::endl;
        return false;
    }
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "----> Error: failed to open the input file " << filename << std::endl;
        return false;
    }

    
std::string firstword, rest, str, type;
while (input >> firstword) {
    
        if (input.eof()) break;

        if(firstword.size() !=0  && firstword[0] == ';'){
            getline(input,rest);
            continue;
        }
//-- simulation (Integrator) block
        if(firstword == AbstractSimulation::GetBaseDefaultReadName()) { // "Integrator_Type"
            input >> str >> type;
            if(type == MC_Simulation::GetDefaultReadName()){
                m_pSimulation = new MC_Simulation(this);
            }
            else if(type == Analysis::GetDefaultReadName()){
                input >> str;
                m_pSimulation = new Analysis(this, str);
            }
            else {
                std::cout<<AbstractSimulation::GetErrorMessage(type)<<std::endl;
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
        else if(firstword == "Box_Centering_F")
        {
            double rate;
            input>>str>>rate;
            m_pSimulation->SetCentering(rate);
            getline(input,rest);
        }
        else if(firstword == "Set_Steps")
        {
            int ini,fi;
            input>>str>>ini>>fi;
            m_pSimulation->UpdateInitialStep(ini);
            m_pSimulation->UpdateFinalStep(fi);
            getline(input,rest);
        }
        else if(firstword == "MinfaceAngle")
        {
            double min_angle;
            input>>str>>min_angle;
            m_pSimulation->SetMinAngle(min_angle);
            getline(input,rest);
        }
        else if(firstword == "Min_Max_Lenghts")
        {
            double min,max;
            input>>str>>min>>max;
            m_pSimulation->SetMinMaxLength(min, max);


            getline(input,rest);
        }
        else if(firstword == "Temperature"){
            
            double beta,delta_beta;
            input>>str>>beta>>delta_beta;
            m_pSimulation->SetBeta(beta,delta_beta);

        }
//---- end of simulation class block
//---- Visualization file
        else if(firstword == AbstractVisualizationFile::GetBaseDefaultReadName()) {
            input>>str>>type;
            if(type == WritevtuFiles::GetDefaultReadName()){  // VTUFileFormat
                double period;
                std::string foldername;
                input >> foldername>> period;
                m_pVisualizationFile = new WritevtuFiles(this, period, foldername);
            }
            getline(input,rest);
        }
//---- End Visualization file
//---- position update
        else if(firstword == AbstractVertexPositionIntegrator::GetBaseDefaultReadName()) {
            input>>str>>type;
            if(type == EvolveVerticesByMetropolisAlgorithm::GetDefaultReadName()){  // MetropolisAlgorithm
                double rate_surf,rate_edge,dr;
                input >> rate_surf >> rate_edge >> dr;
                m_pVertexPositionIntegrator = new EvolveVerticesByMetropolisAlgorithm(this, rate_surf, rate_edge, dr);
            }
            else if(type == EvolveVerticesByKineticMonteCarlo::GetDefaultReadName()){  // KineticMonteCarlo
                double Dv,dt;
                input >>  Dv >> dt;
                m_pVertexPositionIntegrator = new EvolveVerticesByKineticMonteCarlo(this, Dv, dt);
            }
           else if (type == EvolveVerticesByMetropolisAlgorithmWithOpenMPType1::GetDefaultReadName()){  // MetropolisAlgorithmWithOpenMPType1
               double rate_surf,rate_edge,dr;
               input >> rate_surf >> rate_edge >> dr;
#ifdef _OPENMP
               m_pVertexPositionIntegrator = new EvolveVerticesByMetropolisAlgorithmWithOpenMPType1(this, rate_surf, rate_edge, dr);
#else
               std::cout<<"---> warnning: OpenMP was not detected, we will use the serial code "<<type<<"\n";
               m_pVertexPositionIntegrator = new EvolveVerticesByMetropolisAlgorithm(this, rate_surf, rate_edge, dr);
#endif
            }
            else{
                std::cout<<"---> error: unknown method for vertex move "<<type<<"\n";
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
//---- end //
//---- AlexanderMove
        else if(firstword == AbstractAlexanderMove::GetBaseDefaultReadName()) {
                input>>str>>type;
                if(type == AlexanderMoveByMetropolisAlgorithm::GetDefaultReadName()){  // MetropolisAlgorithm
                    double rate;
                    input >> rate;
                    m_pAlexanderMove = new AlexanderMoveByMetropolisAlgorithm(this, rate);
                }
                else if(type == AlexanderMoveByMetropolisAlgorithmWithOpenMP::GetDefaultReadName()){  // MetropolisAlgorithm
                    double rate;
                    input >> rate;
#ifdef _OPENMP
                    m_pAlexanderMove = new AlexanderMoveByMetropolisAlgorithmWithOpenMP(this, rate);
#else
               std::cout<<"---> warnning: OpenMP was not detected, we will use the serial code "<<type<<"\n";
                    m_pAlexanderMove = new AlexanderMoveByMetropolisAlgorithm(this, rate);
#endif
                    
                    
                }
                else{
                    std::cout<<"---> error: unknown method for Alexander move "<<type<<"\n";
                    m_NumberOfErrors++;
                    return false;
                }
            getline(input,rest);
        }
//---- end //
    //---- inclsuion move
            else if(firstword == AbstractInclusionPoseIntegrator::GetBaseDefaultReadName()) {
                    input>>str>>type;
                    if(type == InclusionPoseUpdateByMetropolisAlgorithm::GetDefaultReadName()){  // MetropolisAlgorithm
                        double rate_kawa, rate_angle;
                        input >> rate_kawa>> rate_angle;
                        m_pInclusionPoseIntegrator  = new InclusionPoseUpdateByMetropolisAlgorithm(this, rate_kawa, rate_angle);  // Initialize InclusionPoseIntegrator
                    }
                    else if(type == InclusionPoseUpdateByMetropolisAlgorithmOpenMP::GetDefaultReadName()){  // MetropolisAlgorithm
                        double rate_kawa, rate_angle;
                        input >> rate_kawa>> rate_angle;
#ifdef _OPENMP
                        m_pInclusionPoseIntegrator  = new InclusionPoseUpdateByMetropolisAlgorithmOpenMP(this, rate_kawa, rate_angle);  // Initialize InclusionPoseIntegrator
#else
               std::cout<<"---> warnning: OpenMP was not detected, we will use the serial code "<<type<<"\n";
                        m_pInclusionPoseIntegrator  = new InclusionPoseUpdateByMetropolisAlgorithm(this, rate_kawa, rate_angle);  // Initialize InclusionPoseIntegrator
#endif
                    }
                    else{
                        std::cout<<"---> error: unknown method for Alexander move "<<type<<"\n";
                        m_NumberOfErrors++;
                        return false;
                    }
                getline(input,rest);
            }
    //---- end //
//---- vector field
            else if(firstword == AbstractVectorFieldsRotationMove::GetBaseDefaultReadName()) {
                    input>>str>>type;
                    if(type == VectorFieldsRotationByMetropolisAlgorithm::GetDefaultReadName()){  // MetropolisAlgorithm
                        double rate, dr;
                        input >> rate;
                        getline(input,rest);
                        std::vector<std::string> str_dr = Nfunction::Split(rest);
                        if(str_dr.size() == 0){
                            m_pVectorFieldsRotationIntegrator  = new VectorFieldsRotationByMetropolisAlgorithm (this, rate);  // Initialize InclusionPoseIntegrator
                        }
                        if(str_dr.size() != 0){
                            dr = Nfunction::String_to_Double(str_dr[0]);
                            m_pVectorFieldsRotationIntegrator  = new VectorFieldsRotationByMetropolisAlgorithm (this, rate, dr);  // Initialize InclusionPoseIntegrator
                        }
                    }
                    else{
                        std::cout<<"---> error: unknown method for Alexander move "<<type<<"\n";
                        m_NumberOfErrors++;
                        return false;
                    }
            }
// --------LocalStretching
            else if(firstword == AbstractLocalStretching::GetBaseDefaultReadName()) {
                input>>str>>type;
                if(type == AnisotropicStretchingByNematicField::GetDefaultReadName()){  // MetropolisAlgorithm
                    double kp,kn;
                    input >> kp >> kn;
                    m_pAbstractLocalStretching = new AnisotropicStretchingByNematicField(kp,kn);
                }
            }
//---- Volume_Constraint data
        else if(firstword == AbstractVolumeCoupling::GetBaseDefaultReadName())   { // Volume_Constraint
            input >> str >> type;
            if(type != "No"){
                m_pVAHCalculator->MakeVolumeActive();
                m_pVAHCalculator->MakeAreaActive();
            }
            if(type == VolumeCouplingSecondOrder::GetDefaultReadName()) { // SecondOrderCoupling
                
                double dp,k,vt;
                input>>dp>>k>>vt;
                m_pVolumeCoupling = new VolumeCouplingSecondOrder(m_pVAHCalculator, dp,  k, vt);
            }
            else if(type == Apply_Osmotic_Pressure::GetDefaultReadName()){ // Osmotic_Pressure
                double gamma,P0;
                input>>gamma>>P0;
                m_pVolumeCoupling = new Apply_Osmotic_Pressure(m_pVAHCalculator, gamma, P0 );
            }
            else if(type == "No"){
            }
            else {
                std::cout<<"---> error: unknown volume coupling type "<<type<<"\n";
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
//---- end Volume_Constraint
//-----  start BondedPotentialBetweenVertices
        else if(firstword == AbstractBondedPotentialBetweenVertices::GetBaseDefaultReadName())   { // BondedPotentialBetweenVertices
            input >> str >> type;

            if(type == EmptyBonds::GetDefaultReadName()){
                // it is already set to 
            }
            else if(type == HarmonicBondsList::GetDefaultReadName()) { //
                
                std::string file;
                input>>file;
             m_pBondedPotentialBetweenVertices = new HarmonicBondsList(this, file);
            }
            else {
                std::cout<<AbstractBondedPotentialBetweenVertices::GetErrorMessage(type)<<"\n";
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
    
//-----  start nonBonded
            else if(firstword == AbstractNonbondedInteractionBetweenVertices::GetBaseDefaultReadName())   { // BondedPotentialBetweenVertices
                input >> str >> type;

                if(type == NoNonbondedInteractionBetweenVertices::GetDefaultReadName()){
                    // it is already set to
                    getline(input,rest);
                }
                else if(type == PolarInteractionBetweenEdgesVertices::GetDefaultReadName()) { //
                    
                    getline(input,rest);
                    m_NonbondedInteractionBetweenVertices = new PolarInteractionBetweenEdgesVertices(this, rest);
                }
                else if(type == InteractionBetweenInclusionsIn3D::GetDefaultReadName()) { //
                    
                    getline(input,rest);
                    m_NonbondedInteractionBetweenVertices = new InteractionBetweenInclusionsIn3D(this, rest);
                }
                else {
                    std::cerr<<AbstractNonbondedInteractionBetweenVertices::GetErrorMessage(type)<<"\n";
                    m_NumberOfErrors++;
                    return false;
                }
            }
//---- global curvature
        else if(firstword == AbstractGlobalCurvature::GetBaseDefaultReadName()) {

            input>>str>>type;
            
            if(type != "No"){
                m_pVAHCalculator->MakeGlobalCurvatureActive();
                m_pVAHCalculator->MakeAreaActive();
            }
            
            if(type == CouplingGlobalCurvatureToHarmonicPotential::GetDefaultReadName()) {
                double k,gc0;
                input>>k>>gc0;
                m_pCoupleGlobalCurvature = new CouplingGlobalCurvatureToHarmonicPotential(m_pVAHCalculator,k,gc0);
            }
            else if(type == NoGlobalCurvature::GetDefaultReadName()) {
            }
            else {
                std::cout<<AbstractGlobalCurvature::GetErrorMessage(type)<<"\n";
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
            
        }
// -- global area coupling
        else if(firstword == "TotalAreaCoupling"){
            input>>str>>type;
            if(type != "No"){
                m_pVAHCalculator->MakeAreaActive();
            }
            if(type ==  CouplingTotalAreaToHarmonicPotential::GetDefaultReadName()){
                double gamma , k0;
                input>>k0>>gamma;
                m_pTotalAreaCoupling = new CouplingTotalAreaToHarmonicPotential(m_pVAHCalculator,k0,gamma);
            }
            else if(type == "No") {
            }
            else {
                std::cout<<"---> error: unknown method for global curvature: "<<type<<" \n";
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
// dynamic box
        else if(firstword == AbstractDynamicBox::GetBaseDefaultReadName()){
            int period = 0;
            double force = 0;
            std::string direction;
            input >> str >> type ;

            if (type == PositionRescaleFrameTensionCoupling::GetDefaultReadName()) {
                    
                input >> period >> force >> direction;
                m_pDynamicBox = new PositionRescaleFrameTensionCoupling(period, force, direction, this);
            }
            else if (type == PositionRescaleIsotropicFrameTensionCouplingWithOpenMP::GetDefaultReadName()) {
                input >> period >> force >> direction;
#ifdef _OPENMP
                m_pDynamicBox = new PositionRescaleIsotropicFrameTensionCouplingWithOpenMP(period, force, direction, this);
#else
               std::cout<<"---> warnning: OpenMP was not detected, we will use the serial code "<<type<<"\n";
                m_pDynamicBox = new PositionRescaleFrameTensionCoupling(period, force, direction, this);
#endif
                
            }
            else if (type == PositionRescaleAnisotropicFrameTensionCoupling::GetDefaultReadName()) {
                double force_1 = 0;
                double force_2 = 0;
                double force_3 = 0;
                input >> period >> force_1 >> force_2 >> force_3 >> direction;
                m_pDynamicBox = new PositionRescaleAnisotropicFrameTensionCoupling(period, force_1, force_2, force_3, direction, this);
            }
            else if (type == BoxSizeCouplingToHarmonicPotential::GetDefaultReadName()) {
                double k = 0;
                double a0 = 0;
                input >> period >> k >> a0  >> direction;
                m_pDynamicBox = new BoxSizeCouplingToHarmonicPotential(period, k, a0, direction, this);
            }
            else if (type == BoxSizeCouplingToHarmonicPotentialAnisotropic::GetDefaultReadName()) {
                double kx = 0;
                double ky = 0;
                double kz = 0;
                double a0x = 0;
                double a0y = 0;
                double a0z = 0;

                input >> period >> kx >> ky >> kz >> a0x >> a0y >> a0z >> direction;
                m_pDynamicBox = new BoxSizeCouplingToHarmonicPotentialAnisotropic(period, kx, ky, kz, a0x , a0y, a0z, direction, this);
            }
            else if(type == NoBoxChange::GetDefaultReadName()) {
                    // has been assigned above
            }
            else{
                std::cout<<"--> error: unknown dynamic box method: "<<type<<std::endl;
                m_NumberOfErrors++;
                return false;

            }
            // Consume remaining input line
            getline(input, rest);
        }
// end dynamic box
// open edge treatment
        else if(firstword == AbstractOpenEdgeEvolution::GetBaseDefaultReadName()){ // "OpenEdgeEvolution"
            
            // OpenEdgeEvolution =  EvolutionWithConstantVertex period rate
            input >> str >> type;
            if (type == OpenEdgeEvolutionWithConstantVertex::GetDefaultReadName()) { // "EvolutionWithConstantVertex"
                int period = 0;
                double rate = 0;
                input >>  period>> rate;
                m_pOpenEdgeEvolution = new OpenEdgeEvolutionWithConstantVertex(period, rate, this);
            }
            else if(type == NoEvolution::GetDefaultReadName()){
                
            }
            else {
                
                std::cout<<"---> error: unknown Open Edge Evolution type: "<<type<<"\n";
            }
            getline(input,rest);

        }
// end open edge treatment
        else if(firstword == AbstractDynamicTopology::GetBaseDefaultReadName()){
            input >> str >> type;
            if (type == Three_Edge_Scission::GetDefaultReadName() ) {
                int period = 0;
                input >>  period;

                m_pDynamicTopology = new Three_Edge_Scission(period, this);
            }
            else if(type == ConstantTopology::GetDefaultReadName()){
                m_pDynamicTopology = new ConstantTopology;
            }
            else{
                std::cout<<AbstractDynamicTopology::GetErrorMessage(type);
            }
            // Consume remaining input line
            getline(input, rest);
        }
// InclusionConversion and ActiveTwoStateInclusion
        else if(firstword == AbstractInclusionConversion::GetBaseDefaultReadName()) {
            input >> str >> type;
            
            if(type == ActiveTwoStateInclusion::GetDefaultReadName() ){ // "ActiveTwoStateInclusion"
                
                std::string inctype1, inctype2;
                double ep1,ep2,persentage,gama;
                input>> inctype1>> inctype2>> ep1>> ep2>> persentage>> gama;
                m_pInclusionConversion = new ActiveTwoStateInclusion(ep1, ep2, persentage, gama, inctype1, inctype2);

            }
            else if(type == NoInclusionConversion::GetDefaultReadName()){ // No
                 m_pInclusionConversion = new NoInclusionConversion;
            }
            else{
                 std::cout<<"---> error: unknown Inclusion Conversion method "<<std::endl;
                 m_NumberOfErrors++;
                 return false;
            }
                
            getline(input,rest);
        }
// end InclusionConversion
// boundry condition
        else if(firstword == AbstractBoundary::GetBaseDefaultReadName() ){ // " Boundary "
            input >> str >> type;
            if(type == TwoFlatParallelWall::GetDefaultReadName()){   // " "TwoFlatParallelWall" "
                double thickness;
                char direction;
                input>>thickness>>direction;
                m_pBoundary = new TwoFlatParallelWall(this, thickness,direction);
            }
            else if(type == EllipsoidalShell::GetDefaultReadName()){   // " "EllipsoidalShell" "
                double thickness, r, a, b, c;
                input >> thickness >> r >> a >> b >> c;
                
                m_pBoundary = new EllipsoidalShell(this, thickness, r, a, b, c);
            }
            else if(type == EllipsoidalCore::GetDefaultReadName()){   // " "EllipsoidalShell" "
                double  r, a, b, c;
                input >>  r >> a >> b >> c;
                
                m_pBoundary = new EllipsoidalCore(this, r, a, b, c);
            }
            else {
                std::cout<<"---> error: unknown Boundary type: "<<type<<std::endl;
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
        else if(firstword == "NonequilibriumCommands"){
            getline(input,rest);
            m_pNonequilibriumCommands->LoadCommand(rest);
        }
// end boundry condition
//ConstraintBetweenGroups
        else if(firstword == AbstractApplyConstraintBetweenGroups::GetBaseDefaultReadName()) { // ConstraintBetweenGroups
            
            input >> str >> type;
            if(type == HarmonicPotentialBetweenTwoGroups::GetDefaultReadName()){ // HarmonicPotentialBetweenTwoGroups
                // ConstraintBetweenGroups  = HarmonicPotentialBetweenTwoGroups 10 0.1 2000 Group1 Group2 0 1 1
                double k,l0,nx,ny,nz;
                std::string g1name,g2name;
                input>>k>>l0>>g1name>>g2name>>nx>>ny>>nz;
                m_pApplyConstraintBetweenGroups = new HarmonicPotentialBetweenTwoGroups(this, k, l0, g1name, g2name, nx, ny, nz);
            }
            else if(type == NoConstraint::GetDefaultReadName()){ // No
                
            }
            else{
                std::cout<<" error---> unknown method of Constraint Between Groups "<<std::endl;
            }
            getline(input,rest);

        }
// end ConstraintBetweenGroups
//energy
        else if( firstword == AbstractEnergy::GetBaseDefaultReadName() ) {
            input >> str >> type;
            if(type == Energy::GetDefaultReadName()){ // "FreeDTS1.0_FF"
                m_pEnergyCalculator = new Energy(this);  // Initialize EnergyCalculator
            }else if(type == HMFFEnergy::GetDefaultReadName()){
                std::string data;
                getline(input, data);
                m_pEnergyCalculator = new HMFFEnergy(this, data);
            }
            else {
                std::cout<<" error---> unknown force field type: "<<type<<std::endl;
                m_NumberOfErrors++;
                return false;
            }
        }
//-- end of energy
   //--- extenal force on fields and inclsuions
        else if(firstword == AbstractExternalFieldOnVectorFields::GetBaseDefaultReadName()) {  // ConstantFieldOnVectorFields
            
            input >> str >> type;
            getline(input,rest);
            if(type == ConstantExternalFieldOnVectorFields::GetDefaultReadName() ){
                
                m_pExternalFieldOnVectorFields = new ConstantExternalFieldOnVectorFields(rest);
            }
            else if(type == NoExternalFieldOnVectorFields::GetDefaultReadName() ){
                m_pExternalFieldOnVectorFields = new NoExternalFieldOnVectorFields;

            }
        }
        else if(firstword == AbstractExternalFieldOnInclusions::GetBaseDefaultReadName() ) { // "ConstantField"
            
            input >> str >> type;
            if(type == ConstantExternalField::GetDefaultReadName() ){
                double k,x,y,z;
                input>>k>>x>>y>>z;
                m_pExternalFieldOnInclusions = new ConstantExternalField(k,x,y,z);
            }
            else if(type == ConstantExternalFieldOnOneInclusionType::GetDefaultReadName() ){
                double k,x,y,z;
                std::string inctype;
                input>>inctype>>k>>x>>y>>z;
                m_pExternalFieldOnInclusions = new ConstantExternalFieldOnOneInclusionType(inctype, k,x,y,z);
            }
            else if(type == "No"){
                
            }
            else{
                std::cout<<" unknown External Field On Inclusions method "<<std::endl;
                m_NumberOfErrors++;
                return false;
            }

            getline(input,rest);

        }
        else if(firstword == AbstractVertexAdhesionToSubstrate::GetBaseDefaultReadName() ) { // "ConstantField"
            
            input >> str >> type;
            if(type == SphericalVertexSubstrate::GetDefaultReadName() ){

                getline(input,rest);
                m_pVertexAdhesionToSubstrate = new SphericalVertexSubstrate(rest);
            }
            else if(type == FlatVertexSubstrate::GetDefaultReadName() ){

                getline(input,rest);
                m_pVertexAdhesionToSubstrate = new FlatVertexSubstrate(rest);
            }
            else if(type == "No"){
                getline(input,rest);
            }
            else{
                std::cout<<" unknown AdhesionToSubstrate method "<<std::endl;
                m_NumberOfErrors++;
                getline(input,rest);
                return false;
            }


        }
//-------
        else if( firstword == "MC_Moves" ) {
            
            getline(input, str);
            std::vector<std::string> data = Nfunction::split(str);

            // Check if the correct number of moves is provided
            const int expectedMoves = 5;
            if (data.size() != expectedMoves + 1) { // +1 to account for the first element is just = sign
                std::cout << "---> error: the provided number of MC moves does not cover all the moves. It should be " << expectedMoves << " numbers" << std::endl;
                return false;
            }
            m_pVertexPositionIntegrator->SetMoveRate(Nfunction::String_to_Double(data[1]), Nfunction::String_to_Double(data[2]));
            m_pAlexanderMove->SetMoveRate(Nfunction::String_to_Double(data[3]));
            m_pInclusionPoseIntegrator->SetMoveRate(Nfunction::String_to_Double(data[4]),Nfunction::String_to_Double(data[5]));

        }
        else if(firstword == "FreezingAGroup")
        {
            input>>str>>str;
            m_pVertexPositionIntegrator->UpdateFreezGroupName(str);
            getline(input,rest);
        }
        else if(firstword == AbstractForceonVerticesfromInclusions::GetBaseDefaultReadName()) { //InclusionInducedForceOnVertex
            input >> str >> type;
            if(type == Constant_NematicForceFromAnInclusionType::GetDefaultReadName()){  // Constant_NematicForce
                double f0; // force value
                std::string inc_type;
                input>>f0>>inc_type;
                m_pForceonVerticesfromInclusions = new Constant_NematicForceFromAnInclusionType(f0,inc_type);
            }
            else if(type == Constant_NematicForce::GetDefaultReadName()){  // Constant_NematicForce
                double f0; // force value
                input>>f0;
                m_pForceonVerticesfromInclusions = new Constant_NematicForce(f0);
            }
            else if(firstword == NoForce::GetDefaultReadName()){  // No
                
            }
            else {
                std::cout<<" unknown Force From Inclusions method "<<std::endl;
                m_NumberOfErrors++;
                return false;
            }
            getline(input,rest);
        }
//---- force from vector fields on the vertices
        else if(firstword == AbstractForceonVerticesfromVectorFields::GetBaseDefaultReadName()) {
            input >> str >> type;
            getline(input,rest);
            if(type == Constant_NematicForceByVectorFields::GetDefaultReadName()){  // Constant_NematicForce
                m_pForceonVerticesfromVectorFields = new Constant_NematicForceByVectorFields(rest);
            }
            else if(type == NoVFForce::GetDefaultReadName()){
                
            }
            else {
                std::cout<<" unknown Force From vector fields method "<<std::endl;
                m_NumberOfErrors++;
                return false;
            }
        }
    //---- force from vector fields on the vertices
        else if(firstword == AbstractForceonVertices::GetBaseDefaultReadName()) {
                input >> str >> type;
                getline(input,rest);
                if(type == UserDefinedForceonVertices::GetDefaultReadName()){  // Constant_NematicForce
                    m_pForceonVertices = new UserDefinedForceonVertices(this, rest);
                }
                else if(type == NoForceonVertices::GetDefaultReadName()){
                    
                }
                else {
                    std::cout<<" unknown Force on vertices type has unknown  method "<<std::endl;
                    m_NumberOfErrors++;
                    return false;
                }
        }
        else if(firstword == "TimeSeriesData_Period") { // TimeSeriesData
            int period;
            input>>str>>period;
            m_pTimeSeriesDataOutput->UpdatePeriod(period);
            getline(input,rest);
        }
//-- State class variable
        else if(firstword == "Run_Tag"){
            
                input>>str>>m_GeneralOutputFilename;
                getline(input,rest);
            }
        else if(firstword == "FreezeInclusion")
        {
            std::string inc_freez;
            input>>str>>inc_freez;
            m_pInclusionPoseIntegrator->SetFreezeTypeName(inc_freez);
            getline(input,rest);
        }
        else if(firstword == "Restart_Period")
        {
            int period;
            input>>str>>period;
            m_pRestart->UpdatePeriod(period);
            getline(input,rest);
        }
        else if(firstword == "TopologyFile")
        {
            input>>str>>m_TopologyFile;
            getline(input,rest);
        }
        else if(firstword == "Seed")
        {
            int seed;
            input>>str>>seed;
            m_RandomNumberGenerator = new RNG(seed);
            getline(input,rest);
        }
        else if(firstword == "Kappa")
        {
            double k,kg,c0;
            input>>str>>k>>kg>>c0;
            m_pEnergyCalculator->SetSurfRigidity(k,kg,c0);
            getline(input,rest);
        }
        else if(firstword == "Edge_Parameters")
        {
            double lambda,kg,kn;
            input>>str>>lambda>>kg>>kn;
            m_pEnergyCalculator->SetEdgeRigidity(lambda,kg,kn);
            //=== send it to force field class
            getline(input,rest);
        }
        else if(firstword == "VertexArea")
        {
            double a,b,c,d;
            input>>str>>a>>b>>c>>d;
            
            if(!m_pEnergyCalculator->SetSizeCoupling (a,b,c,d)){
                m_NumberOfErrors++;
                return false;
            }
        }
        else if(firstword == "Voxel_Size"){
            double lx,ly,lz;
            input>>str>>lx>>ly>>lz;
            m_pVoxelization->UpdateVoxelSize(lx, ly,  lz);
            getline(input,rest);
        }
        else if(firstword == AbstractNonbinaryTrajectory::GetBaseDefaultReadName()) {
            
            input>>str>>type;
            if(type == Traj_tsi::GetDefaultReadName()){ // "TSI"
                int period;
                std::string tsiFolder_name;
                input>>tsiFolder_name>>period;
                m_pNonbinaryTrajectory  = new Traj_tsi(this, period, tsiFolder_name);
            }
            getline(input,rest);

        }
        else if(firstword == "ParallelReplica")
        {
            // ParallelReplica = Parallel_Tempering  rate  minbeta    maxbeta
            std::string type;
            input>>str>>type;
            getline(input,rest);

            if(type == "No"|| type == "no"|| type == "NO")
                m_Parallel_Replica.State = false;
            else
                m_Parallel_Replica.State = true;
                        
            m_Parallel_Replica.Type = type;
            m_Parallel_Replica.Data = rest;

        }
        else if(firstword == BTSFile::GetDefaultReadName() ){ // "OutPutTRJ_BTS"
            int periodic, precision;
            std::string filename;
            input>>str>>periodic>>precision>>filename;
            m_pBinaryTrajectory = new BTSFile(this,periodic,precision,filename);
            getline(input,rest);

        }
        else if(firstword == "INCLUSION")
        {
            break;
        }
        else {
            if(firstword.at(0)!=';') {
                std::cout<<"Error: bad keyword in the input file *** "<<firstword<<" ***\n";
                return false;
            }
            getline(input,rest);
        }
    }
    
    // read inclusion
    ReadInclusionType(input);
    
    input.close();
    
    return true;
}
void State::HelpMessage(){
    
    std::cout<<"----------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<"---------------------- weria.pezeshkian@gmail.com ------------------"<<"\n";
    std::cout<<"------------------------------ FreeDTS -----------------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"---------------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    // Print example command
    std::cout << "Example: DTS -f Input.dts -top topol.q -restart res.res" << std::endl;
    std::cout << "=============================================================================" << std::endl;
    
    // Print options
    std::cout <<"     "<<std::setw(20) << std::left << "Option" << std::setw(20) << "Type" << std::setw(20) << "Default" << "Description" << std::endl;
    std::cout <<"----------------------------------------------------------------------------------" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-in" << std::setw(20) << "string" << std::setw(20) << "Input.dts" << "Input file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-top" << std::setw(20) << "string" << std::setw(20) << "topology.top" << "Topology file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-b" << std::setw(20) << "int" << std::setw(20) << "1" << "Initial time step" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-e" << std::setw(20) << "int" << std::setw(20) << "10" << "Final time step" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-seed" << std::setw(20) << "int" << std::setw(20) << "36723" << "Random number seed" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-defout" << std::setw(20) << "string" << std::setw(20) << "dts" << "Output file prefix" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-ndx" << std::setw(20) << "string" << std::setw(20) << "Index.inx" << "Index file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-nt" << std::setw(20) << "int" << std::setw(20) << "1" << "Total number of threads" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-restart" << std::setw(20) << "string" << std::setw(20) << "NO" << "Restart file name" << std::endl;
    std::cout <<"====================================================================================" << std::endl;
    
    return;
}
bool State::Initialize(){
    /*
     * Function: State::Initialize
     * ---------------------------
     * Initializes the state of the simulation by setting up various components such as the mesh, restart files,
     * log files, trajectory files, and data structures. This function performs the following key tasks:
     *
     * 1. Mesh Initialization:
     *    - Creates a MeshBluePrint object using either a restart file or an input topology file.
     *    - Generates the mesh from the mesh blueprint and updates the voxelization box.
     *
     * 2. Restart File Handling:
     *    - Checks if a restart file is provided and attempts to read it.
     *    - If successful, updates various components to the state from the restart file and opens necessary log, data output, and trajectory files.
     *    - If the restart file read fails or is not provided, uses the input topology file to generate the mesh blueprint and opens necessary files for logging and output.
     *
     * 3. Initializing Data Structures and Components:
     *    - Initializes random number generators, curvature calculations, boundary conditions, dynamic topology, inclusion exchange, and various other integrators and couplings.
     *
     * 4. Logging and Visualization:
     *    - Writes the starting state to the log file and creates an initial visualization frame.
     *
     * 5. Energy Calculations:
     *    - Initializes the energy calculator and updates the total energy of the system.
     *
     * 6. Error and Warning Handling:
     *    - Checks for any errors or warnings during initialization and outputs appropriate messages.
     *    - Terminates the program if there are any errors.
     *
     * Returns:
     *    - true if the initialization is successful.
     */
//-----> Get the mesh
        // Create a MeshBluePrint object
        CreateMashBluePrint Create_BluePrint;
        MeshBluePrint mesh_blueprint;
//---->  Check if the restart file name is provided
        bool restartReadSuccess = false;
        if (!m_RestartFileName.empty()) {
            int step;
            double r_vertex;
            double r_box;

            // Attempt to open the restart file
            std::cout << "---> Note: attempting to open the restart file: " << m_RestartFileName << std::endl;

            // Read the restart file and update the State object to that state
            mesh_blueprint = m_pRestart->ReadFromRestart(m_RestartFileName, step, restartReadSuccess, r_vertex, r_box);

            // Check if the restart file was successfully read
            if (restartReadSuccess) {
                std::cout << "---> Note: Restart file was successfully read" << std::endl;
                m_pVertexPositionIntegrator->UpdateDR(r_vertex);
                m_pDynamicBox->UpdateDR(r_box);
                m_pSimulation->UpdateInitialStep(step+1);
                //----- open log file
                if(!m_pTimeSeriesLogInformation->OpenFile(false)){
                    m_NumberOfErrors++;
                }
                if(!m_pTimeSeriesDataOutput->OpenFile(false)){
                    m_NumberOfErrors++;
                }
                //---- open the binary trajectory
                if(!m_pBinaryTrajectory->OpenFile(false, 'w')){
                    m_NumberOfErrors++;
                }
            }
            else{
                // If failed to read restart file, generate MeshBluePrint from input topology file
                std::cout << "---> Warning: Failed to read restart file, will use topology file" << std::endl;
                m_NumberOfWarnings++;
            }
        }
        if(!restartReadSuccess){ // this for also a situation where m_RestartFileName is empty
            mesh_blueprint = Create_BluePrint.MashBluePrintFromInput_Top(m_InputFileName, m_TopologyFile);
            
            //----- open time series files
            if(!m_pTimeSeriesLogInformation->OpenFile(true)){
                m_NumberOfErrors++;
            }
            if(!m_pTimeSeriesDataOutput->OpenFile(true)){
                m_NumberOfErrors++;
            }
            //---- open the binar trajectory
            if(!m_pBinaryTrajectory->OpenFile(true, 'w')){
                m_NumberOfErrors++;
            }
            // Open folder for non-binary trajectory
            if(! m_pNonbinaryTrajectory->OpenFolder()){
                m_NumberOfErrors++;
            }
            // open folder for Visualization
            if(m_pVisualizationFile->GetDerivedDefaultReadName()== "VTUFileFormat"){
                if(!m_pVisualizationFile->OpenFolder()){
                    m_NumberOfErrors++;
                }
            }

        }
    
       // m_RandomNumberGenerator->Initialize();
        // Generate mesh from the mesh blueprint
        m_Mesh.GenerateMesh(mesh_blueprint);
        m_pVoxelization->SetBox(m_pMesh->GetBox());

        // Update group from index file
        m_pMesh->UpdateGroupFromIndexFile(m_IndexFileName);
//============ activate of all inputs and data structures
        m_pRestart->SetRestartFileName();
//----> calaculate curvature for all vertices
        m_pCurvatureCalculations->Initialize();
//----> set some easy access for integrators
        m_pVertexPositionIntegrator->Initialize();
        m_pAlexanderMove->Initialize();
        m_pInclusionPoseIntegrator->Initialize();
        m_pVectorFieldsRotationIntegrator->Initialize();
//----> boundry of the simulations
        m_pBoundary->Initialize();
        m_pForceonVertices->Initialize();
    
//-----> box change
        m_pDynamicBox->Initialize();
    
//----> inclsuion exchange, active inclsuion exchange
    m_pInclusionConversion->Initialize(this);
    m_pDynamicTopology->Initialize();
    m_pOpenEdgeEvolution->Initialize();
    
    m_pAbstractLocalStretching->Initialize();

    
    m_pVAHCalculator->Initialize(this);
    m_pTotalAreaCoupling->Initialize(this);
    m_pVolumeCoupling->Initialize(this);
    m_pCoupleGlobalCurvature->Initialize(this);
   // m_pForceonVerticesfromInclusions this does not have one
    m_pApplyConstraintBetweenGroups->Initialize();
    m_pSimulation->Initialize();
    m_pBondedPotentialBetweenVertices->Initialize();

//--- now that the system is ready for simulation, we first write the State into the log file and make one vis 
    m_pTimeSeriesLogInformation->WriteStartingState();
    m_pVisualizationFile->WriteAFrame(-m_pVisualizationFile->GetPeriod());
//----> energy class
//---> to get interaction energies
    m_pEnergyCalculator->Initialize(m_InputFileName);
//---> to update each vertex and edge energy. Up to now May 2024, edge energy is not zero when both vertices has inclusions
    m_pEnergyCalculator->UpdateTotalEnergy(m_pEnergyCalculator->CalculateAllLocalEnergy());
    
    //-------->
    m_NonbondedInteractionBetweenVertices->Initialize();

    if(m_NumberOfErrors!=0){
        std::cout<<" There were "<<m_NumberOfErrors<<" errors in the input files "<<std::endl;
        exit(0);
    }
    if(m_NumberOfWarnings!=0){
        std::cout<<" There were "<<m_NumberOfWarnings<<" warning in the input files "<<std::endl;
    }
    else{
        std::cout<<" All input files are valid, and the State has been successfully initialized! "<<std::endl;
    }
    return true;
}
bool State::ReadInclusionType(std::ifstream& input) {
    /*
     * @brief Reads and initializes inclusion types from an input file.
     *
     * This function reads the inclusion type definitions from an input file and initializes
     * the necessary data structures in the State object. It reads various properties for
     * each inclusion type and sets up the mesh to use these types.
     *
     * @param input Reference to the input file stream containing inclusion type definitions.
     * @return true if the inclusion types were successfully read and initialized, false otherwise.
     */
    
    std::string firstword, rest, str1, str2, TypeNames;
    int N, TypeID, NoType;
    double Kappa, KappaG, KappaP, KappaL, C0, C0P, C0N;

    // Store inclusion types in a vector
    std::vector<InclusionType> all_InclusionType;

    // Add a default inclusion type
    InclusionType emptyIncType;
    all_InclusionType.push_back(emptyIncType);

    // Read the header line
    input >> str1 >> NoType >> str2;
    getline(input, rest);
    getline(input, rest); // Discard the header line

    // Check if the header line indicates inclusion type definition
    if (str1 == "Define" || str1 == "define" || str1 == "DEFINE") {
        for (int i = 0; i < NoType; i++) {
            std::string inc_data;
            getline(input, inc_data);
            std::vector<std::string> inclusion_str = Nfunction::split(inc_data);
            if (inclusion_str.size() < 9) {
                std::cout << "---> error: bad definition of inclusion types \n";
                m_NumberOfErrors++;
                return false;
            }
            N = Nfunction::String_to_Int(inclusion_str[0]);
            TypeNames = inclusion_str[1];
            Kappa = Nfunction::String_to_Double(inclusion_str[2]);
            KappaG = Nfunction::String_to_Double(inclusion_str[3]);
            KappaP = Nfunction::String_to_Double(inclusion_str[4]);
            KappaL = Nfunction::String_to_Double(inclusion_str[5]);
            C0 = Nfunction::String_to_Double(inclusion_str[6]);
            C0P = Nfunction::String_to_Double(inclusion_str[7]);
            C0N = Nfunction::String_to_Double(inclusion_str[8]);

            // Parse edge data if available
            double lam = 0, ekg = 0, ekn = 0, ecn = 0;
            if (inclusion_str.size() >= 13) {
                lam = Nfunction::String_to_Double(inclusion_str[9]);
                ekg = Nfunction::String_to_Double(inclusion_str[10]);
                ekn = Nfunction::String_to_Double(inclusion_str[11]);
                ecn = Nfunction::String_to_Double(inclusion_str[12]);
            }

            // Create inclusion type and add to vector
            InclusionType incType(TypeNames, i + 1, N, Kappa/2, KappaG, KappaP/2, KappaL/2, C0, C0P, C0N, lam, ekg, ekn, ecn);
            all_InclusionType.push_back(incType);
        }
    }

    // Set inclusion types in the mesh object
    m_Mesh.m_InclusionType = all_InclusionType;

    // Set pointers to inclusion types
    m_Mesh.m_pInclusionType.clear();
    for (size_t i = 0; i < m_Mesh.m_InclusionType.size(); ++i) {
        m_Mesh.m_pInclusionType.push_back(&m_Mesh.m_InclusionType[i]);
    }

    return true;
}
std::string State::CurrentState(){

    std::string state = " Run_Tag = "+ m_GeneralOutputFilename;
    return state;
}

