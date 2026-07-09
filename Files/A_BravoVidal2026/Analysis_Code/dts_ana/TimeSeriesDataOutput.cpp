

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <sstream>
#include "TimeSeriesDataOutput.h"
#include "State.h"
#include "Nfunction.h"
#include "SimDef.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#ifdef MPI_DETECTED
#include <mpi.h>
#endif
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

void TimeSeriesDataOutput::CloseFile(){
    if (m_TimeSeriesFile.is_open()) {
        m_TimeSeriesFile.flush(); 
        m_TimeSeriesFile.close();
     }
}

void TimeSeriesDataOutput::OpenFileWithoutHeader(std::string filename)
{
    m_TimeSeriesFile.open(filename,std::ios_base::app);  
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

    
    m_TimeSeriesFile<<std::fixed;
    m_TimeSeriesFile<<std::setprecision(Precision);
    
//--> write step
    m_TimeSeriesFile<<step<<"   ";

    
//    m_TimeSeriesFile<<step<<"   "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"   ";

//--->write box side length
    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()) {
            m_TimeSeriesFile<<m_pState->GetAnalysisCalculations()->GetEnergy()<<"   ";
        }
        if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()) {
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetArea()<<"   ";
        }
        if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetVolume()<<"   ";
        }
        if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()) {
            m_TimeSeriesFile <<  m_pState->GetAnalysisCalculations()->GetProjectedArea()<< "   ";

        }
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetMeanCurvature() <<"   ";
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetMeanCurvature2() <<"   ";
        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
             m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetGaussianCurvature() <<"   ";
        }

        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            m_TimeSeriesFile <<  m_pState->GetAnalysisCalculations()->GetThickness()<< "   " ;
        }

        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetInclusionNeighbour() << "   ";
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetInclusionEnergy() << "   ";
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetInclusionMeanCurvature() << "   ";
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetInclusionGaussianCurvature() << "   ";
            m_TimeSeriesFile << m_pState->GetAnalysisCalculations()->GetInclusionMeanCurvature2() << "   ";

            if (m_pState->GetAnalysisVariables()->GetNumberInclusionTypes()>1){
                std::vector<std::vector<double>> inclusionneighbours=m_pState->GetAnalysisCalculations()->GetInclusionNeighbourVector();
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    for (int j=i; j<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); j++){
                        m_TimeSeriesFile << inclusionneighbours[i][j] << "   ";
                    }
                   
                }
                std::vector<double> inclusionmeancurvatures=m_pState->GetAnalysisCalculations()->GetInclusionMeanCurvatureVector();
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << inclusionmeancurvatures[i] << "   ";
                }
                std::vector<double> inclusiongaussiancurvatures=m_pState->GetAnalysisCalculations()->GetInclusionGaussianCurvatureVector();
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << inclusiongaussiancurvatures[i] << "   ";
                }
                std::vector<double> inclusionenergies=m_pState->GetAnalysisCalculations()->GetInclusionEnergyVector();
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << inclusionenergies[i] << "   ";
                }
            }
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

    std::string filename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameGeneralAnalysisFile();

    if (!clearfile) {
        // If it's a restart simulation, check if the energy file matches the restart
        m_TimeSeriesFile.open(filename, std::ios_base::app);
        if (!m_TimeSeriesFile.is_open()) {
            std::cerr << "---> error: Unable to open file in append mode: " << filename << std::endl;
            return false;
        }

        // Read the last line of the file
        std::ifstream inputFile(filename);
        if (inputFile.is_open()) {
            std::string lastLine;
            char ch;

            // Move to the end of the file
            inputFile.seekg(0, std::ios_base::end);
            if (inputFile.tellg() == 0) { // Check if the file is empty
                std::cerr << "---> warning: File is empty. No last line to read." << std::endl;
            } else {
                // Backtrack to find the last line
                std::streamoff fileEnd = static_cast<std::streamoff>(inputFile.tellg());
                for (std::streamoff pos = fileEnd - 1; pos >= 0; --pos) {
                    inputFile.seekg(pos, std::ios_base::beg);
                    inputFile.get(ch);
                    if (ch == '\n' && pos != fileEnd - 1) {
                        break; // Stop if a newline is found, and it's not the end
                    }
                    lastLine.insert(lastLine.begin(), ch);
                }

                // Parse the first column of the last line
                std::istringstream lineStream(lastLine);
                int firstColumnValue;
                lineStream >> firstColumnValue;

                if (lineStream.fail()) {
                    std::cerr << "---> warning: Failed to parse the first column from the last line." << std::endl;
                } else {
                    m_pInitialStepRestart=firstColumnValue;
                    // Use firstColumnValue if needed (e.g., for validation)
                }
            }
            inputFile.close();
        } else {
            std::cerr << "---> error: Unable to open file for reading: " << filename << std::endl;
            return false;
        }

    } else {
        // Open the file, clearing it if necessary
        m_TimeSeriesFile.open(filename);
        if (!m_TimeSeriesFile.is_open()) {
            std::cerr << "---> error: Unable to open file: " << filename << std::endl;
            return false;
        }
        // Write the header if it's not a restart
        m_TimeSeriesFile << " ## Frame ";
        if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()) {
            m_TimeSeriesFile << " Energy ";
        }
        if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()) {
            m_TimeSeriesFile << " Area  ";

        }
        if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
            m_TimeSeriesFile << " Volume  ";

        }
        if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()) {
            m_TimeSeriesFile << " Area_p";

        }
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_TimeSeriesFile << " MeanCurvature MeanCurvature2";
            

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_TimeSeriesFile << " GaussianCurvature";
        }

        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            m_TimeSeriesFile << " Thickness";
        }

        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
            m_TimeSeriesFile << " InclusionNeighbour InclusionEnergy InclusionMeanCurvature InclusionGaussianCurvature InclusionMeanCurvature2";

            if (m_pState->GetAnalysisVariables()->GetNumberInclusionTypes()>1){
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    for (int j=i; j<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); j++){
                        m_TimeSeriesFile << " InclusionNeighbourType" + Nfunction::Int_to_String(i) + Nfunction::Int_to_String(j);
                    }
                }

                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << " InclusionMeanCurvatureType" + Nfunction::Int_to_String(i);
                }
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << " InclusionGaussianCurvatureType" + Nfunction::Int_to_String(i);
                }
                for (int i=0; i<m_pState->GetAnalysisVariables()->GetNumberInclusionTypes(); i++){
                    m_TimeSeriesFile << " InclusionEnergyType" + Nfunction::Int_to_String(i);
                }
            }


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
        std::cerr << "Error: Unable to open file, is it here? Why? " << filename << std::endl;
        return false;
    }

    // Create a temporary file for writing
    #ifndef MPI_DETECTED
    std::string NameTemporaryFile="temp.txt";
    #endif

    #ifdef MPI_DETECTED
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string NameTemporaryFile=std::string("temp") + "_" + Nfunction::Int_to_String(rank)+".txt";
    #endif


    std::ofstream tempFile(NameTemporaryFile);
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
    if (std::rename(NameTemporaryFile.c_str(), filename.c_str()) != 0) {
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

void TimeSeriesDataOutput::SetCustomFileName(const std::string& filename) {
    m_customFileName = filename;
}


//Do I need one function that closes the file, and then opens it again?
//IDEA (to be consulted with ChatGPT) --> A function that closes the file, give the file a new name
//and closes the file again... Because, do I have to close it before opening it?
//In that case, I will need to, first of all, make sure that both ranks are closing the file.
//Then, we can open it. Maybe I should set that somehow?