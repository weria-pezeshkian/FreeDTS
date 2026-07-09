#ifndef READTRAJTSIFILES_H
#define READTRAJTSIFILES_H

#include <vector>
#include <string>
#include "SimDef.h"

class ReadTrajTSI {
public:
    ReadTrajTSI(const std::string& folderPath, const std::string& givenStructure);
    bool ValidateFiles();  // Checks the files and populates the vector
    const std::vector<std::string>& GetFilePaths() const;  // Returns the vector of file paths
    int GetNumberOfFrames() const;  // Returns the number of frames
    std::vector<int> GetFrameList();

private:
    std::string m_FolderPath;
    std::string m_FileName;
    std::vector<std::string> m_FilePaths;
    std::vector<int> m_FrameList;  // New member to store the N values
    int m_NumberOfFrames;  // New member to store the number of frames

    bool IsValidFileName(const std::string& fileName, int& index) const;  // Helper function to validate file names
};

#endif // READTRAJTSIFILES_H

