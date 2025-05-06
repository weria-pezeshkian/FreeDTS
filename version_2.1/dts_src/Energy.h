#if !defined(AFX_Energy_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Energy_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Inclusion_Interaction_Map.h"
#include "AbstractEnergy.h"
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
class Energy : public AbstractEnergy {
public:
    Energy(State* pState);
	 ~Energy();

public:
    inline std::string GetDerivedDefaultReadName() {return "FreeDTS1.0_FF";}
    inline static std::string GetDefaultReadName() {return "FreeDTS1.0_FF";}
    double CalculateAllLocalEnergy();   ///
    double SingleVertexEnergy(vertex *p);
    double TwoInclusionsInteractionEnergy(links *);
    double TwoVectorFieldInteractionEnergy(int vf_layer, links * p_edge);
    double CalculateVectorFieldMembraneBindingEnergy(VectorField* p_vf, vertex *p_vertex);
    double CalculateVectorFieldMembraneBindingEnergy(vertex *p_vertex);
private:
    double SurfVertexBendingAndStretchingEnergy(vertex * pver);
    double EdgeVertexBendingAndStretchingEnergy(vertex * pver);
    std::string CurrentState();

    // types of interactions
    double Geo_Theta(vertex *v1, vertex *v2);  // Calculate the angle between two vectors after parallel transport
    double AngleDiff_ParallelTransport(Vec3D &d1, Vec3D &d2, links* p_link);
    double F10(vertex *v1, vertex *v2,std::vector<double>);
    double F2(vertex *v1, vertex *v2,std::vector<double>);
    double F11(vertex *v1, vertex *v2,std::vector<double>);
    double InteractionFunction(double N2, double A, double B, double theta);
    double InteractionFunctionFull(double N2, double A, double B, double theta, links*);
    double InteractionFive(double N2, double A, double B, double C,  links*);
    double Filament_int(double A, double B, double C, vertex* p_v1, vertex* p_v2);

private:
    State* m_pState;
    Vec3D  &m_Box;
    double m_Angle3D;
    double m_Angle2D;




};


#endif
