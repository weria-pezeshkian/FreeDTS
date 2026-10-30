

#include <stdio.h>
#include <fstream>
#include "Inclusion_Interaction_Map.h"
#include "Nfunction.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An object to read pair interaction info of the inclusion type interactions
 This object should be optimised later might even be changed .....
 */
Inclusion_Interaction_Map::Inclusion_Interaction_Map(std::string inputfilename) {

//============================================================
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

    
//---> open the input file
    std::ifstream input(inputfilename.c_str());
    if (!input.is_open()) {
        std::cerr << "Error: Failed to open input file: " << inputfilename << std::endl;
        return;
    }
    std::string line;
    bool hasInteraction = false;
    while (getline(input, line)) {
        if (line == "Inclusion-Inclusion-Int") {
            hasInteraction = true;
            readInteractionData(input);
            break;
        }
    }
    if (!hasInteraction) {
        std::string warning = "Warning: No inclusion interaction energy is mentioned in the input file. They are all set to zero.";
        std::cout << warning << std::endl;
    }
}

Inclusion_Interaction_Map::~Inclusion_Interaction_Map()
{
   
}

bool Inclusion_Interaction_Map::readInteractionData(std::ifstream& input) {
    
    std::string name;
    while (true) {
        int i, j, n;
        double a, b;
        input >> i;
        if (input.eof()) break;
        
        input >> j >> n;
        getline(input, name);
        std::vector<std::string> var = Nfunction::split(name);
        PairInt pint;
        pint.FunctionType = n;
        for (std::vector<std::string>::iterator it = var.begin(); it != var.end(); ++it) {
            if ((*it) == ";") break;
            pint.Varibale.push_back(Nfunction::String_to_Double(*it));
        }

        if (n == 10 && pint.Varibale.size() != 7) {
            std::cerr << "Error: Function type 10 needs 7 parameters but only " << pint.Varibale.size() << " are provided." << std::endl;
            return false; 
        }

        m_MAT[i][j] = pint;
        m_MAT[j][i] = pint;
    }
    
    return true;
}
PairInt Inclusion_Interaction_Map::GetPairInt(int i, int j) {
#if TEST_MODE == Enabled
    if (i >= Inclusion_Type_Number || j >= Inclusion_Type_Number || i < 0 || j < 0) {
        std::cerr << "---> Error: Number of inclusion types exceeds the maximum limit. ";
        std::cerr << "Please adjust the limit in the SimDef.h file." << std::endl;
        std::cerr << "    --: Currently, the limit is set to " << Inclusion_Type_Number << std::endl;
    }
#endif
    return m_MAT[i][j];
}
void Inclusion_Interaction_Map::Print(){
    
    std::cout<<"-----------------------------\n";
    for (int i=0;i<Inclusion_Type_Number;i++) {
        for (int j=0;j<Inclusion_Type_Number;j++) {
            std::vector<double> Varibale  = (m_MAT[i][j]).Varibale;
            for (int n=0;n<Varibale.size();n++)
                std::cout<<Varibale[n]<<"  ";
            std::cout<<"\n";

            
        }
    }

    std::cout<<"-----------------------------\n";

    return;
}
