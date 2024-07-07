#if !defined(AFX_State_H_INCLUDED_)
#define AFX_State_H_INCLUDED_
/*
     State Class: 2024

     Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
     Copyright (c) Weria Pezeshkian

     Description:
     The State class manages the overall state of the simulation. It initializes various components,
     reads input files, explores command-line arguments, and handles the distribution of tasks based on inputs provided.

     The class coordinates the initialization of simulation parameters, such as time steps, input files, and restart options.
     It also sets up integrators, constraints, energy calculators, and a pointer to simulations  for the simulation process.

     The class serves as a central hub for initializing and managing simulation components, ensuring the proper setup
     and execution of the simulation according to the provided inputs.

     Public Methods:
     - State(): Default constructor.
     - State(std::vector<std::string> argument): Constructor that takes command-line arguments.
     - ~State(): Destructor.

     Private Methods:
     - ExploreArguments(std::vector<std::string> &argument): Method to explore and process command-line arguments.
     - ReadInputFile(std::string file): Method to read input files and initialize simulation parameters.
     - HelpMessage(): Method to display a help message with usage instructions.

     Member Variables:
     - Various pointers to components such as TimeSeriesDataOutput, Restart, Traj_tsi, WritevtuFiles, etc.
     - Parameters for simulation setup such as input file names, time steps, thread count, etc.
     - Flags and options for controlling simulation behavior.

     Dependencies:
     - Requires the presence of various supporting classes such as TimeSeriesDataOutput, Restart, Traj_tsi, WritevtuFiles, etc.
     - Relies on external input files for simulation parameters.

     Usage:
     - Construct an instance of the State class to initialize and manage the simulation state.
     - Call appropriate methods to read input files, explore command-line arguments, and execute the simulation.

     Example:
     State stateInstance(arguments);
     stateInstance.ExploreArguments(arguments);
     stateInstance.ReadInputFile(inputFile);
     stateInstance.RunSimulation();

     Notes:
     - Ensure that input files are present and correctly formatted for proper initialization.
     - Command-line arguments should be provided in the correct format to enable proper configuration of the simulation.

 */
#include "SimDef.h"
#include "MESH.h"
#include "Nfunction.h"
//--- system evolution
#include "AbstractVertexPositionIntegrator.h"
#include "AbstractAlexanderMove.h"
#include "AbstractInclusionPoseIntegrator.h"
#include "AbstractVectorFieldsRotationMove.h"
#include "EvolveVerticesByMetropolisAlgorithm.h"
#include "AlexanderMoveByMetropolisAlgorithm.h"
#include "InclusionPoseUpdateByMetropolisAlgorithm.h"
#include "VectorFieldsRotationByMetropolisAlgorithm.h"
//---
//--- I/O
#include "AbstractNonbinaryTrajectory.h"
#include "AbstractVisualizationFile.h"
#include "AbstractBinaryTrajectory.h"
#include "TimeSeriesDataOutput.h"
#include "TimeSeriesLogInformation.h"
#include "Traj_tsi.h"
#include "Restart.h"
#include "WritevtuFiles.h"
#include "BTSFile.h"
//--- energy calcuation method
#include "AbstractEnergy.h"
#include "Energy.h"
//-- curvature
#include "AbstractCurvature.h"
#include "CurvatureByShapeOperatorType1.h"
//-- dynamic box
#include "AbstractDynamicBox.h"
#include "PositionRescaleFrameTensionCoupling.h"
//-- dynamic topology
#include "AbstractDynamicTopology.h"
#include "Three_Edge_Scission.h"
//-- open edge treatment
#include "AbstractOpenEdgeEvolution.h"
#include "OpenEdgeEvolutionWithConstantVertex.h"
//-- volume coupling
#include "AbstractVolumeCoupling.h"
#include "Apply_Osmotic_Pressure.h"
#include "VolumeCouplingSecondOrder.h"
//-- global curvature
#include "AbstractGlobalCurvature.h"
#include "CouplingGlobalCurvatureToHarmonicPotential.h"
//--- total area coupling
#include "AbstractTotalAreaCoupling.h"
#include "CouplingTotalAreaToHarmonicPotential.h"
//--- Constraint Between vertex Groups
#include "AbstractApplyConstraintBetweenGroups.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
//--- Inclusion exchange
#include "AbstractInclusionConversion.h"
#include "ActiveTwoStateInclusion.h"
//--- force from inc to ver
#include "AbstractForceonVerticesfromInclusions.h"
#include "AbstractForceonVerticesfromVectorFields.h"
#include "Constant_NematicForce.h"
#include "Constant_NematicForceByVectorFields.h"
//--- interaction with external fields
#include "AbstractExternalFieldOnVectorFields.h"
#include "AbstractExternalFieldOnInclusions.h"
#include "ConstantExternalField.h"
#include "ConstantExternalFieldOnVectorFields.h"
//--- rigid boundries
#include "AbstractBoundary.h"
#include "RigidWallTypes.h"
#include "Voxelization.h"
//--- simulation
#include "AbstractSimulation.h"
#include "MC_Simulation.h"
//--- accessory objects
#include "RNG.h"
#include "VAHGlobalMeshProperties.h"
#include "InclusionType.h"


struct ParallelReplicaData {  // data structure for turning on and off certain moves
    ParallelReplicaData(){State = false;}
    ParallelReplicaData(bool state){State = state;}
    ~ParallelReplicaData(){}
    bool State;
    std::string Type;
    std::string Data;
};


class State
{
public:
    
	State(std::vector <std::string> argument);
    State();
    ~State();
    
  //  friend class Three_Edge_Scission;

    
//-- standard Integrators
inline AbstractAlexanderMove                *GetAlexanderMove()                 {return m_pAlexanderMove;}
inline AbstractVertexPositionIntegrator     *GetVertexPositionUpdate()                  {return m_pVertexPositionIntegrator;}
inline AbstractInclusionPoseIntegrator      *GetInclusionPoseUpdate()                {return m_pInclusionPoseIntegrator;}
inline AbstractVectorFieldsRotationMove     *GetVectorFieldsRotationUpdate()                {return m_pVectorFieldsRotationIntegrator;}
//---- I/O managment
inline Restart                   *GetRestart()                                      {return m_pRestart;}
inline TimeSeriesDataOutput      *GetTimeSeriesDataOutput()                         {return m_pTimeSeriesDataOutput;}
inline TimeSeriesLogInformation  *GetTimeSeriesLog()                            {return m_pTimeSeriesLogInformation;}
inline AbstractNonbinaryTrajectory   *GetNonbinaryTrajectory()                    {return m_pNonbinaryTrajectory;}
inline AbstractVisualizationFile     *GetVisualization()                           {return m_pVisualizationFile;}
inline AbstractBinaryTrajectory      *GetBinaryTrajectory()                         {return m_pBinaryTrajectory;}

//---- energy and curvature
inline AbstractEnergy               *GetEnergyCalculator()                          {return m_pEnergyCalculator;}
inline AbstractCurvature            *GetCurvatureCalculator()                          {return m_pCurvatureCalculations;}
//----
inline AbstractApplyConstraintBetweenGroups *GetApplyConstraintBetweenGroups()                    {return m_pApplyConstraintBetweenGroups;}

//---- algorithm mangments
inline AbstractTotalAreaCoupling        *GetTotalAreaCoupling()                         {return m_pTotalAreaCoupling;}
inline AbstractVolumeCoupling           *GetVolumeCoupling()                            {return m_pVolumeCoupling;}
inline AbstractGlobalCurvature          *GetGlobalCurvature()                           {return m_pCoupleGlobalCurvature;}
inline AbstractForceonVerticesfromInclusions *GetForceonVerticesfromInclusions()    {return m_pForceonVerticesfromInclusions;}
inline AbstractForceonVerticesfromVectorFields *GetForceonVerticesfromVectorFields()    {return m_pForceonVerticesfromVectorFields;}
inline AbstractExternalFieldOnVectorFields *GetExternalFieldOnVectorFields()        {return m_pExternalFieldOnVectorFields;}
inline AbstractExternalFieldOnInclusions *GetExternalFieldOnInclusions()        {return m_pExternalFieldOnInclusions;}

inline VAHGlobalMeshProperties              *GetVAHGlobalMeshProperties()        {return m_pVAHCalculator;}
//---- supplementary integrators
inline AbstractDynamicBox               *GetDynamicBox()                                {return m_pDynamicBox;}
inline AbstractDynamicTopology          *GetDynamicTopology()                           {return m_pDynamicTopology;}
inline AbstractOpenEdgeEvolution        *GetOpenEdgeEvolution()                         {return m_pOpenEdgeEvolution;}
inline AbstractInclusionConversion      *GetInclusionConversion()                       {return m_pInclusionConversion;}


inline AbstractBoundary                 *GetBoundary()                                  {return m_pBoundary;}
//---- System energy and voxels
inline MESH                     *GetMesh()                                      {return m_pMesh;}  //
inline Voxelization<vertex>     *GetVoxelization()                              {return m_pVoxelization;}
inline AbstractSimulation           *GetSimulation()                                {return m_pSimulation;};
//--- accessory objects
inline  RNG            *GetRandomNumberGenerator()                 const { return m_RandomNumberGenerator; }
//--- some constant variables
inline std::vector <std::string> GetCommandLineArgument()                       {return m_Argument;}
inline std::string               GetRunTag()                                    {return m_GeneralOutputFilename;}
inline int                       GetThreads_Number()                                {return m_Total_no_Threads;}
inline ParallelReplicaData       GetParallelReplicaData()                           {return m_Parallel_Replica;}
std::string CurrentState();

static void HelpMessage();              // writes a help message
bool Initialize(); // makes all the objects ready for simulations, it will open the files ...
    void UpdateRunTag(std::string runtag){
        m_GeneralOutputFilename = runtag;
    }
private:
    bool ReadInputFile(std::string inputfile);    // updates variables based on data in the inputfile
    bool ExploreArguments(std::vector<std::string> &argument);
    bool ReadInclusionType(std::ifstream& input);
private:
    AbstractApplyConstraintBetweenGroups *m_pApplyConstraintBetweenGroups;
    AbstractInclusionConversion* m_pInclusionConversion;
    Restart *m_pRestart;
    AbstractForceonVerticesfromInclusions *m_pForceonVerticesfromInclusions;
    AbstractForceonVerticesfromVectorFields *m_pForceonVerticesfromVectorFields;

    AbstractExternalFieldOnVectorFields *m_pExternalFieldOnVectorFields;
    AbstractExternalFieldOnInclusions *m_pExternalFieldOnInclusions;
//--- Integrators of different degree of freedom
    AbstractAlexanderMove               *m_pAlexanderMove;
    AbstractVertexPositionIntegrator    *m_pVertexPositionIntegrator;
    AbstractInclusionPoseIntegrator     *m_pInclusionPoseIntegrator;
    AbstractVectorFieldsRotationMove    *m_pVectorFieldsRotationIntegrator;
//--- IO file managment
    TimeSeriesDataOutput            *m_pTimeSeriesDataOutput;
    TimeSeriesLogInformation        *m_pTimeSeriesLogInformation;
    AbstractNonbinaryTrajectory         *m_pNonbinaryTrajectory;
    AbstractBinaryTrajectory            *m_pBinaryTrajectory;
    AbstractVisualizationFile           *m_pVisualizationFile;
//---- algorithm mangments
    AbstractDynamicBox            *m_pDynamicBox;
    AbstractDynamicTopology       *m_pDynamicTopology;
    AbstractOpenEdgeEvolution     *m_pOpenEdgeEvolution;
    
    VAHGlobalMeshProperties       *m_pVAHCalculator;
    AbstractVolumeCoupling        *m_pVolumeCoupling;
    AbstractGlobalCurvature       *m_pCoupleGlobalCurvature;
    AbstractTotalAreaCoupling     *m_pTotalAreaCoupling;
    
    /*
     AbstractVolumeCoupling        *m_pVolumeCoupling = new AbstractVolumeCoupling(*m_pVAHCalculator);
     AbstractGlobalCurvature       *m_pCoupleGlobalCurvature = = new AbstractGlobalCurvature(*m_pVAHCalculator);
     AbstractTotalAreaCoupling     *m_pTotalAreaCoupling = new AbstractTotalAreaCoupling (*m_pVAHCalculator);
     
     */
    
    
    AbstractBoundary              *m_pBoundary;
    AbstractCurvature             *m_pCurvatureCalculations;
    AbstractSimulation            *m_pSimulation;
    AbstractEnergy                *m_pEnergyCalculator;

//--- accessory objects
    RNG      *m_RandomNumberGenerator;
//----
    Voxelization<vertex>  *m_pVoxelization;
    MESH                  *m_pMesh;
    MESH                   m_Mesh;

// --- pure members
private:
    std::vector <std::string> m_Argument;
    std::string     m_GeneralOutputFilename; //  a general file flag for specific run
    std::string     m_InputFileName; // name of the topology file, *.top, *.dat *.tsi *.bts
    std::string     m_IndexFileName;            // Name of the index file for group specification
    std::string     m_TopologyFile; // name of the topology file, *.top, *.dat *.tsi *.bts
    std::string     m_RestartFileName; // name of the topology file, *.top, *.dat *.tsi *.bts
    bool m_Targeted_State; // Only relavant for Parallel Tempering by 2022; Which state carries the target temparature
    int m_Total_no_Threads;       // Total no of Threads
    int m_NumberOfErrors;
    int m_NumberOfWarnings;

    
    //--- system variables
        double m_MinVerticesDistanceSquare; //  minimum distance allowed between two vertices  (smaller will results in error)
        double m_MaxLinkLengthSquare;       //  maximum distance allowed between two nighbouring vertices  (larger will results in error)
        double m_MinFaceAngle;              //  minimum angle between the face (smaller will results in error), this is the value of the cos

    ParallelReplicaData m_Parallel_Replica; // an object that includes info about Parallel Tempering method that we are applying
};

#endif
