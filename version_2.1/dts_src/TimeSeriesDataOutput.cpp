

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <sstream>
#include "TimeSeriesDataOutput.h"
#include "State.h"
#include "Nfunction.h"
#include "SimDef.h"

TimeSeriesDataOutput::TimeSeriesDataOutput(){
    m_Periodic = 0;
}
TimeSeriesDataOutput::TimeSeriesDataOutput(State *pState){
    m_pState = pState;
    m_Periodic = 100;
}
TimeSeriesDataOutput::~TimeSeriesDataOutput(){
    if (m_TimeSeriesFile.is_open()) {
        m_TimeSeriesFile.flush(); 
        m_TimeSeriesFile.close();
     }
}
void TimeSeriesDataOutput::UpdatePeriod(int period){
    
    m_Periodic = period;
    return;
}
bool TimeSeriesDataOutput::WriteTimeSeriesDataOutput(int step){
 /*
Writes time series data to the output file.
-param step, The current step of the simulation.
-return True if the data is successfully written, false otherwise.
Note: The function returns true if the data is successfully written,
and false if the periodic condition is not met or if there is an error writing to the file.
*/
    if( m_Periodic == 0 || step%m_Periodic!=0)
        return false;
    
    m_TimeSeriesFile<<std::fixed;
    m_TimeSeriesFile<<std::setprecision(Precision);
    
//--> write step and energy
    m_TimeSeriesFile<<step<<"   "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"   ";

    
//    m_TimeSeriesFile<<step<<"   "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"   ";

//--->write box side length
    if(m_pState->GetDynamicBox()->GetDerivedDefaultReadName()!= NoBoxChange::GetDefaultReadName() ){
        m_TimeSeriesFile<<(*(m_pState->GetMesh()->GetBox()))(0)<<"  ";
        m_TimeSeriesFile<<(*(m_pState->GetMesh()->GetBox()))(1)<<"  ";
        m_TimeSeriesFile<<(*(m_pState->GetMesh()->GetBox()))(2)<<"  ";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->VolumeIsActive()) {
        m_TimeSeriesFile << m_pState->GetVAHGlobalMeshProperties()->GetTotalVolume()<<"  ";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->AreaIsActive()) {
        m_TimeSeriesFile << m_pState->GetVAHGlobalMeshProperties()->GetTotalArea()<<"  ";
    }
    if (m_pState->GetVAHGlobalMeshProperties()->GlobalCurvatureIsActive()) {
        m_TimeSeriesFile << m_pState->GetVAHGlobalMeshProperties()->GetTotalMeanCurvature() <<"  ";
    }
    if (m_pState->GetApplyConstraintBetweenGroups()->GetDerivedDefaultReadName() != NoConstraint::GetDefaultReadName()) {
        m_TimeSeriesFile << m_pState->GetApplyConstraintBetweenGroups()->Output_DataString()<<"  ";
    }
    if (m_pState->GetOpenEdgeEvolution()->GetDerivedDefaultReadName() != "No") {
        m_TimeSeriesFile <<m_pState->GetOpenEdgeEvolution()->GetEdgeSize() <<" ";
    }
    if (m_pState->GetDynamicTopology()->GetDerivedDefaultReadName() != "No") {
        m_TimeSeriesFile <<m_pState->GetDynamicTopology()->GetSurfaceGenus()<<"  ";
    }
    if (m_pState->GetBoundary()->GetDerivedDefaultReadName() != "PBC") {
        m_TimeSeriesFile <<m_pState->GetBoundary()->CurrentStateParameters()<<"  ";
    }
    m_TimeSeriesFile<<std::endl;
    
    return true;
}
bool TimeSeriesDataOutput::OpenFile(bool clearfile) {
    /*
      Opens the time series data output file.
      -param clearfile Flag indicating whether to clear the file before opening.
      -return True if the file is successfully opened, false otherwise.
      This function opens the time series data output file for writing. If the 'clearfile'
      parameter is true, the file is opened in write mode and the header is written. If 'clearfile'
      is false, the file is opened in append mode and checked for compatibility with the restart
      simulation. If the file cannot be opened or there is an error, appropriate error messages
      are printed to stderr.
     */
    std::string filename = m_pState->GetRunTag() + TimeSeriDataExt;

    if (!clearfile) {
        // If it's a restart simulation, check if the energy file matches the restart
        if (!CheckTimeSeriesFile(m_pState->GetSimulation()->GetInitialStep(), filename)) {
            std::cerr << "---> error: energy file does not match with the restart: " << filename << std::endl;
            return false;
        }
        // Append to the file without clearing it
        m_TimeSeriesFile.open(filename, std::ios_base::app);
    } else {
        // Open the file, clearing it if necessary
        m_TimeSeriesFile.open(filename);
        if (!m_TimeSeriesFile.is_open()) {
            std::cerr << "---> error: Unable to open file: " << filename << std::endl;
            return false;
        }
        // Write the header if it's not a restart
        m_TimeSeriesFile << " ## mcstep  energy ";
        if (m_pState->GetDynamicBox()->GetDerivedDefaultReadName() != NoBoxChange::GetDefaultReadName()) {
            m_TimeSeriesFile << " Lx  Ly  Lz ";
        }
        if (m_pState->GetVAHGlobalMeshProperties()->VolumeIsActive()) {
            m_TimeSeriesFile << " Volume  ";

        }
        if (m_pState->GetVAHGlobalMeshProperties()->AreaIsActive()) {
            m_TimeSeriesFile << " Area ";

        }
        if (m_pState->GetVAHGlobalMeshProperties()->GlobalCurvatureIsActive()) {
            m_TimeSeriesFile << " Global_Curvature ";

        }
        if (m_pState->GetApplyConstraintBetweenGroups()->GetDerivedDefaultReadName() != "No") {
            m_TimeSeriesFile << " Constraint_energy ";
        }
        if (m_pState->GetOpenEdgeEvolution()->GetDerivedDefaultReadName() != "No") {
            m_TimeSeriesFile << " edge_size ";
        }
        if (m_pState->GetDynamicTopology()->GetDerivedDefaultReadName() != "No") {
            m_TimeSeriesFile << " surface_genus ";
        }
        if (m_pState->GetBoundary()->GetDerivedDefaultReadName() != "PBC") {
            m_TimeSeriesFile <<" Boundary Info ";
        }
        m_TimeSeriesFile << std::endl;
    }

    return true;
}
bool TimeSeriesDataOutput::FlushFile(){ // the energy file should be flushed first

    m_TimeSeriesFile.flush(); //

    return true;
}
bool TimeSeriesDataOutput::CheckTimeSeriesFile(int ini_step, const std::string& filename) {
    /*
        Checks a time series file for lines with lines written after the last restart has been saved,
        and removes these lines. Somehow update the file.

        Parameters:
            ini_step: The initial step of the new simulations.
            filename: The name of the file to check and modify.

        Returns:
            true if the operation is successful, false otherwise.

        Note:
            This function opens the input file specified by 'filename', reads each line,
            and writes lines with numbers smaller than or equal to 'ini_step' to a temporary file.
            After processing, it replaces the original file with the temporary file.

        If any errors occur during file operations or if an invalid line is encountered
        (i.e., a line without steps), appropriate error messages are printed to standard error.
    */
// Open the input file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return false;
    }

    // Create a temporary file for writing
    std::ofstream tempFile("temp.txt");
    if (!tempFile.is_open()) {
        std::cerr << "Error: Unable to create temporary file" << std::endl;
        inputFile.close();  // Close the input file before returning
        return false;
    }

//----  Read each line from the input file
    std::string line;
    //  read the header
    std::getline(inputFile, line);
    // write the header
    tempFile << line <<std::endl;
    int no_lines = 0;
    while (std::getline(inputFile, line)) {
        no_lines++;
        std::istringstream iss(line);
        int num;
        if (iss >> num) {
            // If the first number on the line is less than or equal to ini_step, write the line to the temp file
            if (num <= ini_step) {
                tempFile << line <<std::endl;  // Use '\n' instead of std::endl for better performance
            }
        } else {
            // Handle invalid lines without a number
            std::cerr << "Warning: Invalid line found in file: " << line << std::endl;
        }
    }
    // Close both input and temporary files
    inputFile.close();
    tempFile.close();

    // Remove the original file
    if (std::remove(filename.c_str()) != 0) {
        std::cerr << "Error: Unable to remove file " << filename << std::endl;
        return false;
    }
    // Rename the temporary file to the original filename
    if (std::rename("temp.txt", filename.c_str()) != 0) {
        std::cerr << "Error: Unable to rename file" << std::endl;
        return false;
    }
    //--- here we check if number of the lines matches the crashed simulation
    if (m_Periodic * (no_lines + 1) < ini_step) {
        
        *(m_pState->GetTimeSeriesLog())<< "---> Warning: An error may have occurred while reading the energy file for restart. \n";
        *(m_pState->GetTimeSeriesLog())<<"       It seems that the input.dts file has been changed compared to the initial file \n";
        return false;
    }

    return true;
}
std::string TimeSeriesDataOutput::CurrentState(){
    
    std::string state = "TimeSeriesData_Period = "+ Nfunction::D2S(m_Periodic);
    return state;
}
