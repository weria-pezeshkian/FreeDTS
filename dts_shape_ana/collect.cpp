#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
#include <map>
#include <iomanip>
#include <valarray>
#include <stdio.h>
double S2D(std::string ConInt);
std::vector<std::string> split(std::string str);
int main(int argc, char* argv[])
{



    std::vector <std::string> argument;
    std::string Temprory;
           for (long i=0;i<argc;i++)
           {
               Temprory.assign(argv[i]);
               argument.push_back(Temprory);
           }
//========================
    std::string infilename = "input.xvg";
    std::string outfilename = "output.xvg";
    int b=0;
    int e=1000000;
    for (long i=1;i<argument.size();i=i+2)
    {
        std::string Arg1 = argument.at(i);
        if(Arg1=="-in")
        {
            infilename = argument.at(i+1);
        }
        else if(Arg1=="-out")
        {
            outfilename = argument.at(i+1);
        }
        else if(Arg1=="-b")
        {
            b = S2D(argument.at(i+1));
        }
        else if(Arg1=="-e")
        {
            e = S2D(argument.at(i+1));
        }
        else
        {
            std::cout<<"error--"<<Arg1<<"\n";
        }
    }
    
    std::ifstream input;
    input.open (infilename.c_str() );
    std::ofstream output;
    output.open (outfilename.c_str() );

    std::string line;
    getline(input,line);
    std::vector<std::vector<double> > Data;
    while(true)
    {
        std::vector<double> temline;
        input>>line;
        int time = S2D(line);
        temline.push_back(time);
        if(input.eof())
            break;
        
        getline(input,line);
        std::vector<std::string> data=split(line);
        
        for (int i=0;i<data.size();i++)
        {
            temline.push_back(S2D(data.at(i)));
        }
        Data.push_back(temline);
    }
    
    std::vector<double>  mean;
    
    for (int i=0;i<(Data.at(0)).size()-1;i++)
        mean.push_back(0);
    
    for (int j=0;j<Data.size();j++)
    for (int i=0;i<(Data.at(0)).size()-1;i++)
    {
        mean.at(i)=mean.at(i)+((Data.at(j)).at(i+1))/double(Data.size());
    }
    
    for (int i=0;i<(Data.at(0)).size()-1;i++)
    {
        std::cout<<mean.at(i)<<"    \n";
    }
    for (int j=0;j<Data.size();j++)
    {
        output<<(Data.at(j)).at(0)<<"   ";
        for (int i=0;i<(Data.at(0)).size()-1;i++)
        {
            double x=(Data.at(j)).at(i+1);
            output<<x<<"   "<<(x-mean.at(i))*(x-mean.at(i))<<"  ";
        }
        output<<"\n";
    }
    

    return 0;
}
double S2D(std::string ConInt)
{
 
     double i =atof( ConInt.c_str() );
    return i;
}
std::vector<std::string> split(std::string str)
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
