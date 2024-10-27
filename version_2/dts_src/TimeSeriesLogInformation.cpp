

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <iostream>
#include "TimeSeriesLogInformation.h"
#include "State.h"
#include "SimDef.h"

TimeSeriesLogInformation::TimeSeriesLogInformation(State *pState){
    m_pState = pState;
}
TimeSeriesLogInformation::~TimeSeriesLogInformation(){
    if (m_TimeSeriesFile.is_open()) {
        m_TimeSeriesFile.flush();
        m_TimeSeriesFile.close();
     }
}
bool TimeSeriesLogInformation::WriteAStringIntoFile(std::string &str){
    
    m_TimeSeriesFile<<str;
    m_TimeSeriesFile<<std::endl;
    
    return true;
}
bool TimeSeriesLogInformation::OpenFile(bool clearfile){
    
    std::string filename = m_pState->GetRunTag();
    filename = filename + ".log";
    
    //-- if it is a restart simulation, then do not clear the file and check for the last line in the file
    if(!clearfile){
        m_TimeSeriesFile.open(filename.c_str(),std::fstream::app);
    }
    else{
        m_TimeSeriesFile.open(filename.c_str());
    }
    
    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "---> error: Unable to open file: " << filename << std::endl;
        return false;
    }
    
    return true;
}
bool TimeSeriesLogInformation::FlushLogFile(){ // the energy file should be flushed first

    m_TimeSeriesFile.flush(); //

    return true;
}
void TimeSeriesLogInformation::WriteStartingState(){
    
    std::vector<std::string> argument = m_pState->GetCommandLineArgument();
    m_TimeSeriesFile<<";--------- this state is initiated by this command ---------------------------------  "<<std::endl;
    for (std::vector<std::string>::iterator it = argument.begin() ; it != argument.end(); ++it){
        m_TimeSeriesFile<<(*it)<<"   ";
    }
    m_TimeSeriesFile<<std::endl;
    m_TimeSeriesFile<<";--------- this part can be used as an input.dts file ---------------------------------  "<<std::endl;
    m_TimeSeriesFile<<m_pState->CurrentState()<<std::endl;
    m_TimeSeriesFile<<";-- abstract classes"<<std::endl;
    m_TimeSeriesFile<<m_pState->GetSimulation()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetBoundary()->CurrentState()<<std::endl;
    
    m_TimeSeriesFile<<m_pState->GetVertexPositionUpdate()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetAlexanderMove()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetInclusionPoseUpdate()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetVectorFieldsRotationUpdate()->CurrentState()<<std::endl;

    m_TimeSeriesFile<<m_pState->GetNonbinaryTrajectory()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetBinaryTrajectory()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetVisualization()->CurrentState()<<std::endl;

    m_TimeSeriesFile<<m_pState->GetCurvatureCalculator()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetEnergyCalculator()->CurrentState()<<std::endl;
    
    m_TimeSeriesFile<<m_pState->GetApplyConstraintBetweenGroups()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetForceonVerticesfromInclusions()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetForceonVerticesfromVectorFields()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetExternalFieldOnVectorFields()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetExternalFieldOnInclusions()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetInclusionConversion()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetVertexAdhesionToSubstrate()->CurrentState()<<std::endl;

    
    m_TimeSeriesFile<<m_pState->GetVolumeCoupling()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetGlobalCurvature()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetTotalAreaCoupling()->CurrentState()<<std::endl;

    m_TimeSeriesFile<<m_pState->GetDynamicBox()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetDynamicTopology()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetOpenEdgeEvolution()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetTimeSeriesDataOutput()->CurrentState()<<std::endl;
    m_TimeSeriesFile<<m_pState->GetRestart()->CurrentState()<<std::endl;

    m_TimeSeriesFile<<";------------------------------------------  "<<std::endl;

}

