#ifndef AFX_TimeSeriesLogInformation_H_INCLUDED_
#define AFX_TimeSeriesLogInformation_H_INCLUDED_

/*
    TimeSeriesLogInformation.h
    
    Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
    Copyright (c) Weria Pezeshkian

    Description:
    This class provides functionality for reading and writing the TimeSeriesLogInformation file.
    It contains methods to write string data into the file, open and flush the log file.
    To write into the file, use the following syntax:
        object << "write something\n";

    Dependencies: State.h
*/

#include <fstream>
#include <string>

class State;

class TimeSeriesLogInformation {
public:
    TimeSeriesLogInformation(State* pstate);
    ~TimeSeriesLogInformation();

public:
    bool OpenFile(bool clear);
    bool FlushLogFile(); // The energy file should be flushed first
    void WriteStartingState();

//----  Overload the << operator to write to the file
    template <typename T>
    TimeSeriesLogInformation& operator<<(const T& data) {
        if (m_TimeSeriesFile.is_open()) {
            m_TimeSeriesFile << data;
        } else {
            std::cerr << "---> error: log file is not open." << std::endl;
        }
        return *this; // Allow chaining of <<
    }
    
    bool WriteAStringIntoFile(std::string &str);

private:
    State *m_pState;
    std::ofstream m_TimeSeriesFile;  // Member variable to hold the file stream
};

#endif
