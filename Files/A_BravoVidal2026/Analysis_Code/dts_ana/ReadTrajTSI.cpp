#include "ReadTrajTSI.h"
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <iostream>  // Add this to use std::cout
#include <algorithm>  // For std::sort


namespace fs = std::filesystem;

ReadTrajTSI::ReadTrajTSI(const std::string& folderPath, const std::string& givenStructure)
    : m_FolderPath(folderPath), m_FileName(givenStructure), m_NumberOfFrames(0) {}

bool ReadTrajTSI::ValidateFiles() {
    m_FilePaths.clear();
    m_FrameList.clear();
    
    if (!fs::exists(m_FolderPath) || !fs::is_directory(m_FolderPath)) {
        return false;  // Folder does not exist or is not a directory
    }

    // First pass: Collect N values from valid filenames
    for (const auto& entry : fs::directory_iterator(m_FolderPath)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().filename().string();
            int index = -1;
            if (IsValidFileName(fileName, index)) {
                
                
                m_FrameList.push_back(index);
            
            }
        }
    }



    if (m_FrameList.empty()) {
        return false;  // No valid files found
    }
   
    // Sort m_FrameList from low to high
    std::sort(m_FrameList.begin(), m_FrameList.end());

    // Set m_NumberOfFrames and resize m_FilePaths accordingly
    m_NumberOfFrames = static_cast<int>(m_FrameList.size());
    m_FilePaths.resize(m_NumberOfFrames);

    // Second pass: Populate m_FilePaths based on m_FrameList
    for (int i = 0; i < m_NumberOfFrames; ++i) {
        int frameIndex = m_FrameList[i];
        std::stringstream ss;
        ss << m_FileName << frameIndex << ".tsi";
        m_FilePaths[i] = m_FolderPath + "/" + ss.str();  // Construct file path
    }

    return true;  // Successfully processed files
}

const std::vector<std::string>& ReadTrajTSI::GetFilePaths() const {
    return m_FilePaths;
}

int ReadTrajTSI::GetNumberOfFrames() const {
    return m_NumberOfFrames;
}

std::vector<int> ReadTrajTSI::GetFrameList() {
    return m_FrameList;
}

bool ReadTrajTSI::IsValidFileName(const std::string& fileName, int& index) const {
    std::string prefix = m_FileName;  // Expect the fileName to start with prefix + "_"
    std::string suffix = ".tsi";

    // Ensure the filename has both the prefix and suffix, and that they are correctly positioned
    if (fileName.size() <= prefix.size() + suffix.size()) {
        return false;
    }

    // Check if the filename starts with the correct prefix and ends with the correct suffix
    if (fileName.compare(0, prefix.size(), prefix) != 0 || fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) != 0) {
        return false;
    }

    // Extract the numeric part between the prefix and suffix
    std::string indexPart = fileName.substr(prefix.size(), fileName.size() - prefix.size() - suffix.size());

    // Ensure that the index part consists only of digits
    if (indexPart.empty() || !std::all_of(indexPart.begin(), indexPart.end(), ::isdigit)) {
        return false;
    }

    // Convert the numeric part to an integer
    try {
        index = std::stoi(indexPart);
    } catch (const std::exception&) {
        return false;  // In case the conversion fails
    }

    return true;
}