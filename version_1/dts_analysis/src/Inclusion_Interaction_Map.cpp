

#include <stdio.h>
#include "Inclusion_Interaction_Map.h"
#include "Nfunction.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An object to read pair interaction info of the inclusion type interactions
 This object should be optimised later might even be changed .....
 */
Inclusion_Interaction_Map::Inclusion_Interaction_Map()
{
}
Inclusion_Interaction_Map::Inclusion_Interaction_Map(std::string inputfilename)
{
for (int i=0;i<Inclusion_Type_Number;i++)
{
    std::vector<PairInt> T;
    PairInt TP;
    TP.FunctionType = 0;
    for (int j=0;j<Inclusion_Type_Number;j++)
    {
        T.push_back(TP);
    }
    m_MAT.push_back(T);
}
#if TEST_MODE == Enabled
    std::cout<<"---->Note: Make sure these interactions are defined at the end of the file -- "<<std::endl;
#endif
    
// Read parameters from the inputfile
    Nfunction f;
    std::string line;
    std::string FileName = inputfilename;
    std::ifstream input;
    input.open(FileName.c_str());
    std::string name;
    bool hasit=false;
    while (true)
    {
        getline(input,name);
        if(input.eof())
        {
            if(hasit==false)
            {
#if TEST_MODE == Enabled
            std::string sms="----> Note: no inclsuion interaction energy is mentioned in the inputfile, they are all set to zero";
            std::cout<<sms<<"\n";
#endif
            }
            break;
        }
        
        if(name=="Inclusion-Inclusion-Int")
        {
            hasit = true;
            while (true)
            {
                int i,j,n;
                double a,b;
                input>>i;
                if(input.eof())
                    break;
                else
                {
                    input>>j>>n;
                    getline(input,name);
                    std::vector<std::string> var  = f.split(name);
                    std::vector<double> doublevar;
                    PairInt pint;
                    for (std::vector<std::string>::iterator it = var.begin() ; it != var.end(); ++it)
                    {
                        if((*it)==";")
                            break;
                        else
                            doublevar.push_back(f.String_to_Double(*it));
                        
                    }
                    pint.FunctionType = n;
                    pint.Varibale = doublevar;
                    if(n==10 && doublevar.size()!=7)
                    {
                        std::cout<<"Error: Function type 10 need 7 parameters but "<<doublevar.size()<<" are provided \n";
                    }
                    (m_MAT.at(i)).at(j) = pint;
                    (m_MAT.at(j)).at(i) = pint;

                }
                
            
            }
        }
        
        
    }

}

Inclusion_Interaction_Map::~Inclusion_Interaction_Map()
{
   
}
//=== The efficiency is needed only for this function and how do we access its members. 
PairInt Inclusion_Interaction_Map::GetPairInt (int i,int j)
{
#if TEST_MODE == Enabled
    if(i>Inclusion_Type_Number || j>Inclusion_Type_Number || i<0 || j<0)
    {
        std::string sms="Error: Bead types in the interaction energy";
        std::cout<<sms<<"\n";
    }
#endif
    PairInt Int = (m_MAT.at(i)).at(j);
    return Int;
}

