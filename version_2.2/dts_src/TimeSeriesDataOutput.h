#ifndef AFX_TimeSeriesDataOutput_H_INCLUDED_
#define AFX_TimeSeriesDataOutput_H_INCLUDED_

/*
    TimeSeriesDataOutput.h
    
    Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
    Copyright (c) Weria Pezeshkian

    Description:
    This class provides functionality for reading and writing the TimeSeriesDataOutput file.
    It contains methods to write time series data, open and flush the output file, and check compatibility
    of the file with a restart simulation.
    
    Dependencies: State.h, Nfunction.h, SimDef.h
*/

#include <fstream>

class State;
class TimeSeriesDataOutput {
public:
    TimeSeriesDataOutput();
    TimeSeriesDataOutput(State* pstate);
    ~TimeSeriesDataOutput();
    
    inline int GetPeriodic() const { return m_Periodic; }

    inline static std::string GetDefaultReadName() {return "TimeSeriesData";}

public:
    bool WriteTimeSeriesDataOutput(int step);
    bool OpenFile(bool clearfile);
    bool FlushFile(); // The energy file should be flushed first
    void UpdatePeriod(int period);
    std::string CurrentState();

private:
    State *m_pState;
    int m_Periodic;
    std::ofstream m_TimeSeriesFile;  // Member variable to hold the file stream
    bool CheckTimeSeriesFile(int endstep, const std::string &filename); // This checks if the energy file matches the restart end step
};

#endif
