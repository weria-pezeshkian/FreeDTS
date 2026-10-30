#ifndef AFX_NFUNCTION_H_7F4A21B8_B13C_11D3_BF19_004095086186__INCLUDED_
#define AFX_NFUNCTION_H_7F4A21B8_B13C_11D3_BF19_004095086186__INCLUDED_

#include <string>
#include <vector>
class Vec3D;
class Nfunction {
public:
    // Conversion functions
    static std::string Int_to_String(double value);
    static std::string D2S(double value);
    static int String_to_Int(const std::string& str);
    static double String_to_Double(const std::string& str);
    // String manipulation functions
    static std::string SubstringFromRight(const std::string& input, const char& chr);
    // File operations
    static bool FileExist(const std::string& name);
    // Utility functions
    static bool isEven(int x);
    static bool CopyFile(std::string file1, std::string file2);
    
    static std::vector<std::string> split(const std::string& str);
    static std::vector<std::string> Split(const std::string& str); // new version of it

    static bool OpenFolder(const std::string &foldername);

    static bool CopyBinaryFile(const std::string& file1, const std::string& file2, const std::streamsize bufferSize);
    static void HelpMessage();
    static std::string ConvertSecond2Time(double seconds);
    static bool isValidDoubleNumber(double d);
    
    static double SquarePBCDistanceOfTwoPoint(const Vec3D &P1, const Vec3D &P2, const Vec3D &Box);
    static void ConsolePrint_Note(const std::string &text);
    static void ConsolePrint_Error(const std::string &text);
    static void ConsolePrint_Warning(const std::string &text);

};

#endif

