
#include <fstream>
#include <time.h>
#include <sstream>
#include "MESH.h"
#include "CreateMashBluePrint.h"
#include "State.h"
#include "Traj_tsi.h"
#include "WritevtuFiles.h"
#include "BTSFile.h"
#include "./Registry/FactoryVolumeCoupling.h"
#include "./Registry/FactoryDynamicTopologyMethod.h"
#include "./Registry/FactoryVertexPositionIntegrator.h"
#include "./Registry/FactoryDynamicBox.h"

// new version, all in one. Since 2026
#include "./Registry/FactoryForRegistration.h"

State::State(){
    
}
State::State(std::vector<std::string> argument) :

      m_pMesh(&m_Mesh), // still until the opening of the constractor, this will not be affective
      m_NumberOfErrors(0),
      m_NumberOfWarnings(0),
      m_CanSimulationCall(true),
      m_CanEnergyCall(true),
      m_CanpVPositionIntegratorCall(true),
      m_CanIPoseIntegratorCall(true),
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
      //---- Initialize supplementary integrators
      m_pDynamicBox(new NoBoxChange),  // Initialize DynamicBox
      m_pDynamicTopology(new ConstantTopology),  // Initialize DynamicTopology
      m_pOpenEdgeEvolution(new NoEvolution),  // Initialize OpenEdgeEvolution
      m_pInclusionConversion(new NoInclusionConversion),  // Initialize InclusionConversion

      //--- Initialize accessory objects
      m_pCurvatureCalculations(new CurvatureByShapeOperatorType1(this)),  // Initialize CurvatureCalculations
      m_RandomNumberGenerator(new RNG(1234)),  // Initialize RandomNumberGenerator
      //m_pEnergyCalculator(new Energy(this)),  // Initialize EnergyCalculator

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
    //m_pMesh = &m_Mesh;
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
    std::ifstream input_read(filename);
    if (!input_read) {
        std::cerr << "----> Error: failed to open the input file " << filename << std::endl;
        return false;
    }
std::string firstword, rest, str, type, line;
    
while (std::getline(input_read, line)) {
    
        if (line.empty() || line[0] == ';') continue;
    
        std::istringstream input(line);
        input >> firstword;
        
        if (firstword == AbstractSimulation::GetBaseDefaultReadName() && m_CanSimulationCall ) {
                 if (!RegisterUsingInputFile<FactorySimulationScheme>(input, "SC442878", m_pSimulation)) {
                     return false;
                 }
        }
        else if (firstword == AbstractEnergy::GetBaseDefaultReadName() && m_CanEnergyCall ) {
            if (!RegisterUsingInputFile<FactoryEnergyFunctions>(input, "EF326834", m_pEnergyCalculator)) {
                    return false;
            }
        }
    // open edge treatment
        else if (firstword == AbstractOpenEdgeEvolution::GetBaseDefaultReadName()) {
            if (!RegisterUsingInputFile<FactoryOpenEdgeEvolutionMethod>(input, "OE292838", m_pOpenEdgeEvolution)) {
                      return false;
                  }
        }
    // end open edge treatment
    //---- position update

        else if (firstword == AbstractVertexPositionIntegrator::GetBaseDefaultReadName() && m_CanpVPositionIntegratorCall) {
            if (!RegisterUsingInputFile<FactoryVertexPositionIntegrator>(input, "VP292838", m_pVertexPositionIntegrator)) {
                    return false;
            }
        }
    //---- end position update //
   // dynamic box,
        else if (firstword == AbstractDynamicBox::GetBaseDefaultReadName()) {
            if (!RegisterUsingInputFile<FactoryDynamicBox>(input, "DB392838", m_pDynamicBox)) {
                    return false;
            }
        }
    // end dynamic box
    //---- AlexanderMove
        else if (firstword == AbstractAlexanderMove::GetBaseDefaultReadName()) {
            if (!RegisterUsingInputFile<FactoryAlexanderMove>(input, "AM586838", m_pAlexanderMove)) {
                    return false;
            }
        }
    //---- end //
        else if (firstword == AbstractDynamicTopology::GetBaseDefaultReadName()) {
                if (!RegisterUsingInputFile<FactoryDynamicTopologyMethod>(input, "DT986858", m_pDynamicTopology)) {
                    return false;
                }
        }
    //---- Visualization file
        else if (firstword == AbstractVisualizationFile::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryVisualizationFile>(input, "VF446878", m_pVisualizationFile)) {
                    return false;
            }
        }
    //---- End Visualization file
    //--- ConstraintBetweenGroups
        else if (firstword == AbstractApplyConstraintBetweenGroups::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryConstraintBetweenGroups>(input, "CG346878", m_pApplyConstraintBetweenGroups)) {
                    return false;
            }
        }
    //End
    //---- inclsuion move
        else if (firstword == AbstractInclusionPoseIntegrator::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryInclusionPoseIntegrator>(input, "IP676878", m_pInclusionPoseIntegrator)) {
                    return false;
            }
        }
    // curvature method
        else if (firstword == AbstractCurvature::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryCurvatureMethod>(input, "CM096878", m_pCurvatureCalculations)) {
                    return false;
            }
        }
    // external field on incs
        else if (firstword == AbstractExternalFieldOnInclusions::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryExternalFieldOnInclusions>(input, "EF876378", m_pExternalFieldOnInclusions)) {
                    return false;
            }
        }
        else if (firstword == AbstractVectorFieldsRotationMove::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryVectorFieldsRotationMove>(input, "VR887378", m_pVectorFieldsRotationIntegrator)) {
                    return false;
            }
        }
        else if (firstword == AbstractInclusionConversion::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryInclusionConversionMethod>(input, "IC287378", m_pInclusionConversion)) {
                    return false;
            }
        }
        else if (firstword == AbstractNonbondedInteractionBetweenVertices::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryNonbondedInteractionBetweenVertices>(input, "NB387378", m_NonbondedInteractionBetweenVertices)) {
                    return false;
            }
        }
        else if (firstword == AbstractForceonVerticesfromInclusions::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryForceonVerticesfromInclusions>(input, "FI387378", m_pForceonVerticesfromInclusions)) {
                    return false;
            }
        }
        else if (firstword == AbstractForceonVertices::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryForceonVertices>(input, "FV457378", m_pForceonVertices)) {
                    return false;
            }
        }
        else if (firstword == AbstractExternalFieldOnVectorFields::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryExternalFieldOnVectorFields>(input, "EV37378", m_pExternalFieldOnVectorFields)) {
                    return false;
            }
        }
        else if (firstword == AbstractForceonVerticesfromVectorFields::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryForceonVerticesfromVectorFields>(input, "VF99378", m_pForceonVerticesfromVectorFields)) {
                    return false;
            }
        }
        else if (firstword == AbstractNonbinaryTrajectory::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryNonbinaryTrajectory>(input, "NT00378", m_pNonbinaryTrajectory)) {
                    return false;
            }
        }
        else if (firstword == AbstractBinaryTrajectory::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryBinaryTrajectory>(input, "BT00278", m_pBinaryTrajectory)) {
                    return false;
            }
        }
        else if (firstword == AbstractBoundary::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryBoxBoundary>(input, "BC00291", m_pBoundary)) {
                    return false;
            }
        }
        else if (firstword == AbstractVertexAdhesionToSubstrate::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryVertexAdhesionToSubstrate>(input, "BC00291", m_pVertexAdhesionToSubstrate)) {
                    return false;
            }
        }
        else if (firstword == AbstractBondedPotentialBetweenVertices::GetBaseDefaultReadName()) {
             if (!RegisterUsingInputFile<FactoryBondedPotentialBetweenVertices>(input, "BC00291", m_pBondedPotentialBetweenVertices)) {
                    return false;
            }
        }
    // -- global area coupling
            else if(firstword == AbstractTotalAreaCoupling::GetBaseDefaultReadName()){
                    
                    input >> str >> type;
                    m_pVAHCalculator->MakeAreaActive();
                    m_pTotalAreaCoupling = FactroyTotalAreaCoupling::Instance().Create(type, input, m_pVAHCalculator);

                if(!m_pTotalAreaCoupling) {
                    std::cout << AbstractTotalAreaCoupling::GetRegistryError(type)<< std::endl;
                    m_NumberOfErrors++;
                    return false;
                }
            }
                else if(firstword == AbstractVolumeCoupling::GetBaseDefaultReadName()) {
                    
                    input >> str >> type;
                    
                        m_pVAHCalculator->MakeVolumeActive();
                        m_pVAHCalculator->MakeAreaActive();
                        m_pVolumeCoupling = FactoryVolumeCoupling::Instance().Create(type, input, m_pVAHCalculator);
                        
                        if(!m_pVolumeCoupling) {
                            std::cout << AbstractVolumeCoupling::GetRegistryError(type)<< std::endl;
                            m_NumberOfErrors++;
                            return false;
                        }
                }
        else if(firstword == AbstractGlobalCurvature::GetBaseDefaultReadName()) {

                input>>str>>type;
                m_pVAHCalculator->MakeGlobalCurvatureActive();
                m_pVAHCalculator->MakeAreaActive();

                m_pCoupleGlobalCurvature = FactroyGlobalCurvature::Instance().Create(type, input, m_pVAHCalculator);

            if(!m_pCoupleGlobalCurvature) {
                std::cout<<AbstractGlobalCurvature::GetErrorMessage(type)<<"\n";
                m_NumberOfErrors++;
                return false;
            }
        }
    //---- end //
        else if(HandleSimulationCommands(firstword, input)){
            m_CanSimulationCall = false;
        }
        else if(HandleEnergyCommands(firstword, input)){
            m_CanEnergyCall = false;
        }
        else if(HandleVPositionIntegratorCommands(firstword, input)){
            m_CanpVPositionIntegratorCall = false;
        }
        else if(HandleIPoseIntegratorCommands(firstword, input)){
            m_CanIPoseIntegratorCall = false;
        }
    // this is different, this is loading
        else if(firstword == "NonequilibriumCommands"){
            getline(input,rest);
            m_pNonequilibriumCommands->LoadCommand(rest);
        }


        else if(firstword == "TimeSeriesData_Period") { // TimeSeriesData
            int period;
            input>>str>>period;
            m_pTimeSeriesDataOutput->UpdatePeriod(period);
        }
//-- State class variable
        else if(firstword == "Run_Tag"){
            
                input>>str>>m_GeneralOutputFilename;
            }
        else if(firstword == "Restart_Period")
        {
            int period;
            input>>str>>period;
            m_pRestart->UpdatePeriod(period);
        }
        else if(firstword == "TopologyFile")
        {
            input>>str>>m_TopologyFile;
        }
        else if(firstword == "Seed")
        {
            int seed;
            input>>str>>seed;
            m_RandomNumberGenerator = new RNG(seed);
        }
        else if(firstword == "Voxel_Size"){
            double lx,ly,lz;
            input>>str>>lx>>ly>>lz;
            m_pVoxelization->UpdateVoxelSize(lx, ly,  lz);
            m_pVoxelization->SetDefaultVoxelSize(lx, ly,  lz);
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
        else if(firstword == "INCLUSION")
        {
            break;
        }
        else
        {
            if (!firstword.empty() && firstword[0] != ';')
            {
                std::cerr
                    << "**********************************************************\n"
                    << "-->Error: unknown or unexpected keyword in input file: '"
                    << firstword << "'\n"
                    << "This indicates that the parser does not recognize this command.\n"
                    << "Possible causes:\n"
                    << "  - typo or misspelling in the keyword\n"
                    << "  - unsupported or deprecated parameter name\n"
                    << "  - mismatch between input file version and code\n"
                    << "  e.g., Integrator_Type if declared should come first \n"
                    << "  and., EnergyMethod if declared should come second \n"
                << "**********************************************************";


                return false;
            }
            std::getline(input, rest);
        }
    }
    
    // read inclusion
    ReadInclusionType(input_read);
    
    input_read.close();
    
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
       
        if(!(m_Mesh.GenerateMesh(mesh_blueprint))){
            std::cout<<"\033[1;31m [ERROR]-->\033[0m An error has been found during mesh generation  \n";
            m_NumberOfErrors++;
        }
        else{

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
            m_pInclusionConversion->Initialize();
            m_pDynamicTopology->Initialize();
            m_pOpenEdgeEvolution->Initialize();
    
    
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
    }
    if(m_NumberOfErrors != 0){
        std::cout<<" There were "<<m_NumberOfErrors<<" errors in the input files "<<std::endl;
        return false;
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
bool State::HandleSimulationCommands(const std::string& firstword, std::istream& input){
    
    std::string str;
    if(firstword == "Box_Centering_F")
    {
        double rate;
        input>>str>>rate;
        m_pSimulation->SetCentering(rate);
        return true;
    }
    else if(firstword == "Set_Steps")
    {
        int ini,fi;
        input>>str>>ini>>fi;
        m_pSimulation->UpdateInitialStep(ini);
        m_pSimulation->UpdateFinalStep(fi);
        return true;
    }
    else if(firstword == "MinfaceAngle")
    {
        double min_angle;
        input>>str>>min_angle;
        m_pSimulation->SetMinAngle(min_angle);
        return true;
    }
    else if(firstword == "Min_Max_Lenghts")
    {
        double min,max;
        input>>str>>min>>max;
        m_pSimulation->SetMinMaxLength(min, max);
        return true;

    }
    else if(firstword == "Temperature"){
        double beta,delta_beta;
        input>>str>>beta>>delta_beta;
        m_pSimulation->SetBeta(beta,delta_beta);
        return true;
    }
    else if(firstword == "Check_Mesh_Frequency")
    {
        double rate;
        input>>str>>rate;
        m_pSimulation->CheckMeshFrequency(rate);
        return true;
    }
    
    return false;
}
bool State::HandleEnergyCommands(const std::string& firstword, std::istream& input){
    
    std::string str;
    if(firstword == "Kappa")
    {
        double k,kg,c0;
        input>>str>>k>>kg>>c0;
        m_pEnergyCalculator->SetSurfRigidity(k,kg,c0);
        return true;
    }
    else if(firstword == "Edge_Parameters")
    {
        double lambda,kg,kn;
        input>>str>>lambda>>kg>>kn;
        m_pEnergyCalculator->SetEdgeRigidity(lambda,kg,kn);
        return true;
    }
    else if(firstword == "VertexArea")
    {
        double a,b,c,d;
        input>>str>>a>>b>>c>>d;
        m_pEnergyCalculator->SetSizeCoupling (a,b,c,d);
        return true;
    }
    
    return false;
}
bool State::HandleVPositionIntegratorCommands(const std::string& firstword, std::istream& input){
    
    std::string str;
    if(firstword == "FreezingAGroup")
    {
        input>>str>>str;
        m_pVertexPositionIntegrator->UpdateFreezGroupName(str);
        return true;

    }
    return false;
}
bool State::HandleIPoseIntegratorCommands(const std::string& firstword, std::istream& input){
    
    if(firstword == "FreezeInclusion")
    {
        std::string inc_freez, str;
        input>>str>>inc_freez;
        m_pInclusionPoseIntegrator->SetFreezeTypeName(inc_freez);
        return true;
    }
    
    return false;
}
std::string State::CurrentState(){

    std::string state = " Run_Tag = "+ m_GeneralOutputFilename;
    return state;
}

