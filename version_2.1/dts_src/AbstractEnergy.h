#if !defined(AFX_AbstractEnergy_H)
#define AFX_AbstractEnergy_H
#include <iostream>
#include "Inclusion_Interaction_Map.h"
// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for system energy. Only includes the individual vertex and inclsuion energy
========================================================
*/
class  AbstractEnergy {
public:
    AbstractEnergy(){
        m_Lambda = 0;        // line tension
        m_Kappa_Geo = 0;    // geodesic curvature rigidity for edge V
        m_Kappa_Norm = 0;   // normal curvature rigidty for edge V
        m_kappa = 0;        // bending rigidity
        m_kappa_G = 0;      // gaussian curvature rigidity
        m_SCurvature0 = 0;
        m_Ka = 0;          // coupling for mean area
        m_Area0 = 0;     // mean area of v
        m_Kl =0;        // couple for mean edge length
        m_l0 =0;        // mean edge length
        
    }
    virtual ~ AbstractEnergy(){
        
        delete m_pInt;
    }
    //---- virtual functions
    
    virtual inline double CalculateAllLocalEnergy() = 0 ;
    virtual double TwoInclusionsInteractionEnergy(links *) = 0;
    virtual double TwoVectorFieldInteractionEnergy(int layer, links *) = 0;
    virtual double SingleVertexEnergy(vertex *p) = 0;
    virtual double CalculateVectorFieldMembraneBindingEnergy(VectorField* p_vf, vertex *p_vertex) = 0;
    virtual double CalculateVectorFieldMembraneBindingEnergy(vertex *p_vertex) = 0;
    
    virtual std::string CurrentState() = 0;

    
    void print(){m_pInt->Print();}
    //--- 
    
    
    void Initialize(std::string inputfilename){
        
        m_pInt = new Inclusion_Interaction_Map(inputfilename);
        // Check if m_pInt was successfully allocated
        if (m_pInt == NULL) {
            std::cerr << "Error: Failed to allocate memory for Inclusion_Interaction_Map \n";
        }
        
        return;
    }
   // virtual  double CalaculateEnergyofSingleVertex(vertex * pvertex) = 0;


    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "EnergyMethod";}
    
//---> real functions
    inline double GetEnergy()               const                 {return m_TotalEnergy;}
    
    void UpdateTotalEnergy(double en){
       
        m_TotalEnergy = en;
        return;
    }
    void AddToTotalEnergy(double d_en){
       
        m_TotalEnergy += d_en;
        return;
    }
    
    void SetSurfRigidity(double kappa, double kappa_g, double c0){
        m_kappa = kappa/2;
        m_kappa_G = kappa_g;
        m_SCurvature0 = c0;
        return;
    }
    void SetEdgeRigidity(double lambda, double kappa_geo, double kappa_norm){
        m_Lambda = lambda;
        m_Kappa_Geo = kappa_geo;
        m_Kappa_Norm = kappa_norm;
        return;
    }
    bool SetSizeCoupling(double ka, double a0, double kl, double l0){
        
        if(a0<0 || a0>1 || l0<0 || l0>1){
            std::cout<<"---> error in constant vertex area; gamma is bad; it should be a double number between 0-1 \n";
            return false;
        }
        m_Ka = ka;
        m_Area0 = (1+2*a0)*sqrt(3)/2.0;
        m_Kl = kl;
        m_l0 = sqrt(1+2*l0);
        return true;
    }
private:
    double m_TotalEnergy;
    
protected:
    double m_Lambda;    // line tension
    double m_Kappa_Geo;  // geodesic curvature rigidity for edge V
    double m_Kappa_Norm;  // normal curvature rigidty for edge V
    double m_kappa; // bending rigidity
    double m_kappa_G; // gaussian curvature rigidity
    double m_SCurvature0; // gaussian curvature rigidity
    double m_Ka;          // coupling for mean area
    double m_Area0;     // mean area of v
    double m_Kl;        // couple for mean edge length
    double m_l0;        // mean edge length
    
    Inclusion_Interaction_Map * m_pInt;

};

#endif
