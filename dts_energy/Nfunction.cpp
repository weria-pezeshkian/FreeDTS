#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>
#include <vector>
#include "Nfunction.h" 
Nfunction::Nfunction()
{


}
Nfunction::~Nfunction()
{

}
std::string Nfunction::Int_to_String(double ConInt)
{
    std::string StrCon;  
    std::ostringstream S_mediate;  
    S_mediate << ConInt;  
    StrCon = S_mediate.str();
    return StrCon;
}
void Nfunction::CleanFiles()
{
    
    std::ofstream logfile;
    logfile.open("sim_run.log");
    logfile.close();
}


int Nfunction::String_to_Int(std::string ConInt)
{
 
     int i =atoi( ConInt.c_str() );
    return i;
}


double Nfunction::String_to_Double(std::string ConInt)
{
 
     double i =atof( ConInt.c_str() );
    return i;
}



double Nfunction::Word_Searcher(std::string Word_Name,std::string File_Name)
{
std::ifstream FILE_NAME;
FILE_NAME.open (File_Name.c_str() ,std::ios::binary);

L=0;
while (L==0)
{
FILE_NAME>>WORD;
L=bool((WORD==Word_Name)||FILE_NAME.eof());
}

if (FILE_NAME.eof())
{
std::cout<<"There is no command"<<Word_Name;
return 0;
}
else 
{
FILE_NAME>>NUMBER;
return NUMBER;
}
FILE_NAME.close();
}
//
std::string Nfunction::Word_Searcher_String(std::string Word_Name,std::string File_Name)
{
std::string Front_Word;
std::ifstream FILE_NAME;
FILE_NAME.open (File_Name.c_str() ,std::ios::binary);

L=0;
while (L==0)
{
FILE_NAME>>WORD;
L=bool((WORD==Word_Name)||FILE_NAME.eof());
}

if (FILE_NAME.eof())
{
std::cout<<"this command was not found "<<Word_Name;
return 0;
}
else 
{
FILE_NAME>>Front_Word;
return Front_Word;
}
FILE_NAME.close();
}

int Nfunction::AFactor(int R)
{
        int a=1;
    if(R==0)
    {
            int a=1;
    }

    for (int i=2;i<R+1;i++)
    {
        a=a*i;
    }
    return a;

}

bool Nfunction::isEven(int x)
{
 if (x % 2 == 0)
{
return true;
}
else
{
return false;
}

}
void Nfunction::Write_One_LogMessage(std::string A)
{
    std::ofstream logfile;
    logfile.open("sim_run.log",std::fstream::app);
    logfile<<A<<"\n";
    logfile.close();
   // 
    
}
std::vector<std::string> Nfunction::split(std::string str)
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
float Nfunction::Quicksqrt(float num)
{
    long i;
    float x2, y;
    
    const float threehalfs = 1.5F;
    
    x2 = num * 0.5F;
    y = num;
    
    i= * (long * ) &y;
    i = 0x5f3759df - (i>>1);
    y = * (float *) &i;
    y = y * (threehalfs - (x2*y*y));
    y = y * (threehalfs - (x2*y*y));
    y = y * (threehalfs - (x2*y*y));
    
    return y;
}
