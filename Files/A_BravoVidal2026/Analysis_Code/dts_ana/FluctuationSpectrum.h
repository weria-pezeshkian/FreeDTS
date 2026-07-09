#if !defined(AFX_FluctuationSpectrum_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_FluctuationSpectrum_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
#include <complex>
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"
#include <fstream>
using Complex = std::complex<double>;

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
struct FourierResult;
class State;
class FluctuationSpectrum {
public:
    FluctuationSpectrum(State* pState, int Nx, int Ny,int type);
	 ~FluctuationSpectrum();

public:
    void CalculateSpectrum();
    void OpenOutputStreams(bool clearfile);
    void CloseOutputStreams();
    void CalculateSpectrumIndividual();
    void CalculateSpectrumIndividualComplex();
    void CalculateSpectrumSphericalComplex();
private:
    void GenerateVectorOrder();
    void GenerateVectorOrderSpherical();
    void GenerateZeroAndNonZeroVectorOrder();
    void GenerateZeroAndNonZeroVectorOrderIndividual();
    double FourierTransformNoInclusion(std::vector<std::vector<double>> qvector);
    std::vector<double> FourierTransformNoInclusionIndividual(std::vector<std::vector<double>> qvector);
    std::vector<double> FourierTransformInclusion(std::vector<std::vector<double>> qvector);
    std::vector<std::vector<double>> FourierTransformInclusionIndividual(std::vector<std::vector<double>> qvector);
    void AverageHeight();
    void GenerateMatrixOrder();
    void CenterOfMass();
    FourierResult FourierTransformInclusionIndividualComplex(std::vector<std::vector<double>> qvector);
    std::vector<Complex> FourierTransformNoInclusionIndividualComplex(std::vector<std::vector<double>> qvector);
    FourierResult FourierTransformInclusionSphericalComplex(std::vector<std::vector<int>> qvector);
    void SphericalHarmonic(int l, int m, double theta, double phi, double &real, double &imag);
    double Factorial(int n);
    double AssociatedLegendre(int l, int m, double x);


    


private:
    State* m_pState;
    int m_Nx,m_Ny;
    std::vector<std::vector<int>> m_VectorOrder;
    int m_SpectrumSize;
    std::vector<std::vector<std::vector<int>>> m_MatrixOrder;
    std::vector<std::vector<int>> m_MatrixOrderIndividual;
    double m_AverageHeight;
    double m_AverageInclusionDensity;
    std::ofstream m_QVector,m_HQVector,m_HPQVector,m_PQVector;
    int m_type;
    std::vector<double> m_CenterOfMass;
    double m_radius0;
    int m_lmax;
    double m_Lx,m_Ly,m_Lz;
    
    
    


};


#endif
