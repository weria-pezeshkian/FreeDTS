#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>
#include <vector>
#include <cctype>
#include <iomanip>
#include <string>
#include "SimDef.h"
//#include "Nfunction.h"

std::string Nfunction::Int_to_String(double ConInt) {
    std::ostringstream S_mediate;
    S_mediate << ConInt;
    return S_mediate.str();
}
std::string Nfunction::D2S(double ConInt) {
    std::ostringstream S_mediate;
    S_mediate << ConInt;
    return S_mediate.str();
}
int Nfunction::String_to_Int(const std::string& ConInt) {
    return std::atoi(ConInt.c_str());
}
bool Nfunction::OpenFolder(const std::string &foldername) {
    // Use std::experimental::filesystem to create the directory
    const int dir_err = system(("mkdir -p "+foldername).c_str());
    if (-1 == dir_err){
        std::cout<<"---> ---> Error: Failed to create directory  "<<foldername<<"\n";
        return false;
    }
    return true;
}
double Nfunction::String_to_Double(const std::string& ConInt) {
    return std::atof(ConInt.c_str());
}
bool Nfunction::isEven(int x) {
    return x % 2 == 0;
}
std::vector<std::string> Nfunction::split(const std::string& str)
{
    
    std::vector<std::string> Line;
    
    std::string word;
    bool flag=false;
    for (int i=0;i<str.size();i++)
    {
        if (str.at(i) ==' '||str.at(i) =='\f' || str.at(i) =='\n' || str.at(i) =='\r' || str.at(i) =='\t' || str.at(i) =='\v' )
        {
            if(flag == true)
            {
                Line.push_back(word);
                word.clear();
                flag=false;
            }
        }
        else
        {
            flag = true;
            word.push_back(str.at(i));
        }
    }
    if(flag == true)
    {
        Line.push_back(word);
        word.clear();
        flag=false;
    }
    return Line;
}
bool Nfunction::FileExist (const std::string& name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
std::string Nfunction::SubstringFromRight(const std::string& input, const char& chr) {
    size_t pos = input.find_last_of(chr);
    if (pos != std::string::npos) {
        return input.substr(pos + 1);
    } else {
        // If there is no '/' in the string, return the entire string
        return input;
    }
}
bool Nfunction::CopyFile(std::string file1, std::string file2)
{
    //==============================================================================
    //=========== This function copies file1 into file2
    std::ifstream in(file1.c_str());
    std::ofstream out(file2.c_str());
    std::string str;
    if(in.is_open() && out.is_open())
    {
        while(true)
        {
            getline(in,str);
            if(in.eof())
                break;
            out<<str<<"\n";
        }
    }
    
    //Close both files
    in.close();
    out.close();
    
    return true;
}
std::vector<std::string> Nfunction::Split(const std::string& str)
{
    std::vector<std::string> words; // Vector to store individual words
    std::string word; // Temporary string to store each word
    bool inWord = false; // Flag to indicate if currently reading a word

    // Iterate through each character in the input string
    for (size_t i = 0; i < str.size(); ++i) {
        char c = str[i];

        // Check if the character is a whitespace
        if (std::isspace(c)) {
            // If we were reading a word, add it to the vector
            if (inWord) {
                words.push_back(word);
                word.clear(); // Clear the temporary string
                inWord = false; // Reset the flag
            }
        } else {
            // Add the character to the temporary word string
            word.push_back(c);
            inWord = true; // Set the flag to indicate that we are reading a word
        }
    }

    // If there's a word left after processing the string, add it to the vector
    if (!word.empty()) {
        words.push_back(word);
    }

    return words; // Return the vector of words
}
bool Nfunction::CopyBinaryFile(const std::string& inputFilename, const std::string& outputFilename, const std::streamsize bufferSize) {
    std::ifstream inFile(inputFilename, std::ios::binary);
    std::ofstream outFile(outputFilename, std::ios::binary);

    if (!inFile.is_open()) {
        std::cerr << "Error: Unable to open input file: " << inputFilename << std::endl;
        return false;
    }

    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open output file: " << outputFilename << std::endl;
        return false;
    }
    char buffer[bufferSize];

    while (true) {
        inFile.read(buffer, bufferSize);
        std::streamsize bytesRead = inFile.gcount();

        if (bytesRead > 0) {
            outFile.write(buffer, bytesRead);
        }

        if (inFile.eof()) {
            break;
        }

        if (!inFile) {
            std::cerr << "Error: Failed to read from input file: " << inputFilename << std::endl;
            return false;
        }

        if (!outFile) {
            std::cerr << "Error: Failed to write to output file: " << outputFilename << std::endl;
            return false;
        }
    }

    return true;
}
void Nfunction::HelpMessage() {
    // Print header
    std::cout << "=============================================================================" << std::endl;
    std::cout << "    Welcome to FreeDTS - Dynamically Triangulated Surface Simulation            " << std::endl;
    std::cout << "=============================================================================" << std::endl;
    
    // Print software version
    std::cout << "Version: " << SoftWareVersion << std::endl;
    std::cout << "=============================================================================" << std::endl;
    
    // Print example command
    std::cout << "Example: DTS -f Input.dts -top topol.q -restart res.res" << std::endl;
    std::cout << "=============================================================================" << std::endl;
    
    // Print options
    std::cout << std::left << std::setw(12) << "Option" << std::setw(15) << "Type" << std::setw(20) << "Default" << "Description" << std::endl;
    std::cout << "----------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(12) << "-in" << std::setw(15) << "string" << std::setw(20) << "Input.dts" << "Input file name" << std::endl;
    std::cout << std::left << std::setw(12) << "-top" << std::setw(15) << "string" << std::setw(20) << "topology.top" << "Topology file name" << std::endl;
    std::cout << std::left << std::setw(12) << "-b" << std::setw(15) << "int" << std::setw(20) << "1" << "Initial time step" << std::endl;
    std::cout << std::left << std::setw(12) << "-e" << std::setw(15) << "int" << std::setw(20) << "10" << "Final time step" << std::endl;
    std::cout << std::left << std::setw(12) << "-seed" << std::setw(15) << "int" << std::setw(20) << "36723" << "Random number seed" << std::endl;
    std::cout << std::left << std::setw(12) << "-defout" << std::setw(15) << "string" << std::setw(20) << "dts" << "Output file prefix" << std::endl;
    std::cout << std::left << std::setw(12) << "-ndx" << std::setw(15) << "string" << std::setw(20) << "Index.inx" << "Index file name" << std::endl;
    std::cout << std::left << std::setw(12) << "-nt" << std::setw(15) << "int" << std::setw(20) << "1" << "Total number of threads" << std::endl;
    std::cout << std::left << std::setw(12) << "-restart" << std::setw(15) << "string" << std::setw(20) << "NO" << "Restart file name" << std::endl;
    std::cout << "=============================================================================" << std::endl;
}
std::string Nfunction::ConvertSecond2Time(double Seconds) {
    int seconds = Seconds;
    if(seconds != 0)
    Seconds = Seconds/double(seconds);
    const int daysPerMonth = 30; // Approximate average days per month
    const int hoursPerDay = 24;
    const int minutesPerHour = 60;
    const int secondsPerMinute = 60;

    // Convert seconds to months, days, hours, minutes, and seconds
    int months = seconds / (daysPerMonth * hoursPerDay * minutesPerHour * secondsPerMinute);
    seconds %= (daysPerMonth * hoursPerDay * minutesPerHour * secondsPerMinute);

    int days = seconds / (hoursPerDay * minutesPerHour * secondsPerMinute);
    seconds %= (hoursPerDay * minutesPerHour * secondsPerMinute);

    int hours = seconds / (minutesPerHour * secondsPerMinute);
    seconds %= (minutesPerHour * secondsPerMinute);

    int minutes = seconds / secondsPerMinute;
    seconds %= secondsPerMinute;

    // Create a stringstream to build the output string
    std::stringstream ss;

    // Append non-zero units to the stringstream
    bool printedSomething = false;

    if (months > 0) {
        ss << months << " months";
        printedSomething = true;
    }
    if (days > 0) {
        if (printedSomething) ss << ", ";
        ss << days << " days";
        printedSomething = true;
    }
    if (hours > 0) {
        if (printedSomething) ss << ", ";
        ss << hours << " hours";
        printedSomething = true;
    }
    if (minutes > 0) {
        if (printedSomething) ss << ", ";
        ss << minutes << " minutes";
        printedSomething = true;
    }
    if (seconds > 0 || !printedSomething) {
        if (printedSomething) ss << ", ";
        ss << seconds + Seconds << " seconds";
    }

    // Return the built string
    return ss.str();
}
