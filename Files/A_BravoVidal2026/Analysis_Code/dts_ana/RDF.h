#if !defined(AFX_RDF_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_RDF_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"
#include <fstream>

/*
 * @brief Energy calculation based on FreeDTS1.1 force field.
 *
 * This class is responsible for calculating the energy of a system using the FreeDTS1.1 force field.
 * It includes methods to compute various energy contributions such as single vertex energy, energy
 * due to interactions between two inclusions.

 *
 * @note This class inherits from the AbstractEnergy interface, providing a common interface for energy
 * calculation modules in the simulation framework.
 *
 * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * @copyright Weria Pezeshkian
 */
class State;
class RDF {
public:
    RDF(State* pState,double dr);
	 ~RDF();

public:
    void CalculateRDF();
    void OpenOutputStreams(bool clearfile);
    void CloseOutputStreams();
    void create_distance_vector(double r_min, double r_max, double dr);
    void create_normalization_vector();
    
private:
    State* m_pState;
    std::ofstream m_InclusionCluster;
    int m_pNumberofInclusions;
    double m_Lx;
    double m_Ly;
    double m_dr;
    int m_pNumberofBins;
    std::vector<double> m_Rvector;
    std::vector<double> m_Nvector;
};


#endif
