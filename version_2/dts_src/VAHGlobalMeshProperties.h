#if !defined(AFX_VAHGLOBALMESHPROPERTIES_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_VAHGLOBALMESHPROPERTIES_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_

/*
 * Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 *
 * This class is created in version 1.2 to centralize global variables like volume, area, and total curvature.
 * It provides methods to calculate and retrieve these global properties of a mesh.
 */

#include "SimDef.h"

class vertex;
class triangle;
class links;
class State;
class AbstractVolumeCoupling;
class AbstractTotalAreaCoupling;
class AbstractGlobalCurvature;

class VAHGlobalMeshProperties {
public:
    VAHGlobalMeshProperties();
    virtual ~VAHGlobalMeshProperties();

    // Friend classes to allow them to access private members of this class
    friend class AbstractVolumeCoupling;
    friend class AbstractTotalAreaCoupling;
    friend class AbstractGlobalCurvature;

    // Inline getters for the global properties
    inline double GetTotalVolume() const { return m_TotalVolume; }
    inline double GetTotalArea() const { return m_TotalArea; }
    inline double GetTotalMeanCurvature() const { return m_TotalCurvature; }
    inline bool GetCalculateVAH() const { return m_CalculatedGlobalVariable; }

    // Public methods for calculating contributions to global variables
    void CalculateAVertexRingContributionToGlobalVariables(vertex* p_vertex, double& vol, double& area, double& curvature);
    void CalculateALinkTrianglesContributionToGlobalVariables(links* p_link, double& vol, double& area, double& curvature);
    void CalculateBoxRescalingContributionToGlobalVariables(double lx, double ly, double lz, double& vol, double& area, double& curvature);
    void Initialize(State* pState);

    // Method to update the calculation status of global variables
    void UpdateCalculatedGlobalVariable() {
        m_CalculatedGlobalVariable = true;
    }

    // Method to calculate the volume of a single triangle
    double CalculateSingleTriangleVolume(triangle* pTriangle);

private:
    // Private member variables to store global properties
    double m_TotalVolume;
    double m_TotalArea;
    double m_TotalCurvature; // Delta A = h * m_TotalCurvature = h * Sum [2H_v * A_v]
    bool m_CalculatedGlobalVariable;
    State* m_pState; // Pointer to the state object
};

#endif // AFX_VAHGLOBALMESHPROPERTIES_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
