#ifndef AFX_Inclusion_Interaction_Map_H_9P4B21B8_C13C_5648_BF23_124095086237__INCLUDED_
#define AFX_Inclusion_Interaction_Map_H_9P4B21B8_C13C_5648_BF23_124095086237__INCLUDED_

#include "SimDef.h"
#include "Vec3D.h"

/*
 * @brief Inclusion Interaction Map object
 *
 * This class handles reading and storing pair interaction information
 * for inclusion type interactions.
 *
 * Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */
struct PairInt {
    int FunctionType;
    std::vector<double> Varibale;
};

class Inclusion_Interaction_Map {
public:
    // Constructor to initialize the interaction map from a file
    Inclusion_Interaction_Map(std::string inputfilename);

    // Destructor
    ~Inclusion_Interaction_Map();

    // Get the pair interaction for given indices i and j
    PairInt GetPairInt(int i, int j);
    void Print();

private:
    // Matrix to store pair interactions
    std::vector<std::vector<PairInt> > m_MAT;

    // Function to read one interaction data
    bool readInteractionData(std::ifstream& input);
};

#endif
