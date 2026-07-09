#ifndef ANALYSISCALCULATIONS_H
#define ANALYSISCALCULATIONS_H

#include <vector>
#include <string>
#include "SimDef.h"
#include "Nfunction.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"

class State;
class AnalysisCalculations {
public:

    AnalysisCalculations(State *pState);

    
    double GetArea() {return m_pArea;}

    double GetVolume() {return m_pVolume;}

    
    double GetEnergy() {return m_pEnergy;}

    
    double GetProjectedArea() {return m_pProjectedArea;}

    
    double GetMeanCurvature() {return m_pMeanCurvature;}
    double GetMeanCurvature2() {return m_pMeanCurvature2;}

    
    double GetThickness() {return m_pThickness;}

    
    double GetGaussianCurvature() {return m_pGaussianCurvature;}

    double GetInclusionNeighbour() {return m_pInclusionNeighbour;}
    double GetInclusionEnergy() {return m_pInclusionEnergy;}
    double GetInclusionMeanCurvature() {return m_pInclusionMeanCurvature;}
    double GetInclusionMeanCurvature2() {return m_pInclusionMeanCurvature2;}
    double GetInclusionGaussianCurvature() {return m_pInclusionGaussianCurvature;}
    std::vector<std::vector<double>> GetInclusionNeighbourVector() {return m_pInclusionNeighbourVector;}

    std::vector<double> GetInclusionMeanCurvatureVector() {return m_pInclusionMeanCurvatureVector;}
    std::vector<double> GetInclusionGaussianCurvatureVector() {return m_pInclusionGaussianCurvatureVector;}
    std::vector<double> GetInclusionEnergyVector() {return m_pInclusionEnergyVector;}
    

    void Calculate();

private:

    /*void CalculateArea();
    void CalculateEnergy();
    void CalculateProjectedArea();
    void CalculateMeanCurvature();
    void CalculateThickness();
    void CalculateGaussianCurvature();*/
    void InitializeMemberVariables();
    double CalculateSingleTriangleVolume(triangle *pTriangle);
    void CalculateInclusion();
    void CalculateInclusionNeighbours();


private:

    double m_pGaussianCurvature;
    double m_pProjectedArea;
    double m_pEnergy;
    double m_pArea;
    double m_pThickness;
    double m_pMeanCurvature;
    double m_pMeanCurvature2;
    double m_pVolume;
    double m_pInclusionNeighbour;
    double m_pInclusionEnergy;
    double m_pInclusionMeanCurvature;
    double m_pInclusionMeanCurvature2;
    double m_pInclusionGaussianCurvature;
    std::vector<std::vector<double>> m_pInclusionNeighbourVector;
    std::vector<double> m_pInclusionMeanCurvatureVector;
    std::vector<double> m_pInclusionGaussianCurvatureVector;
    std::vector<double> m_pInclusionEnergyVector;
    std::vector<int> m_pNumberOfInclusionPerTypesVector;
    State *m_pState;
};

#endif // READTRAJTSIFILES_H