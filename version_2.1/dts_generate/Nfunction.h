#if !defined(AFX_NFUNCTION_H_7F4A21B8_B13C_11D3_BF19_004095086186__INCLUDED_)
#define AFX_NFUNCTION_H_7F4A21B8_B13C_11D3_BF19_004095086186__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Some functions ...
 */
class Nfunction
{
public:
Nfunction();
~Nfunction();
    
    
    public:
    int AFactor(int R);
    std::string Int_to_String(double ); //   A function for making a string from a intger
    int String_to_Int(std::string);
    double String_to_Double(std::string);
    double Word_Searcher(std::string, std::string); // A function get a word and a filename and looking for word into the file
    std::string Word_Searcher_String(std::string, std::string); // A function get a word and a filename and looking for word into the file
    void Write_One_LogMessage(std::string);
    void Write_One_ErrorMessage(std::string);
    void CleanFiles();
    std::vector<std::string> split(std::string str);
    bool FileExist (const std::string& name);
    bool isEven(int x);
    private :
    int L;
    double NUMBER;
    std::string WORD;
    float Quicksqrt(float num);


};





#endif
