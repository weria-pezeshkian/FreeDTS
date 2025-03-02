#if !defined(AFX_AbstractBondedPotentialBetweenVertices_H)
#define AFX_AbstractBondedPotentialBetweenVertices_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for BondedPotentialBetweenVertices
========================================================
*/
class  AbstractBondedPotentialBetweenVertices {
public:
    AbstractBondedPotentialBetweenVertices(){
        
    }
    virtual ~ AbstractBondedPotentialBetweenVertices(){
        
    }
    virtual void Initialize() = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
 //   virtual double GetBondEnergyOfVertex(vertex *pvertex) = 0;
    virtual double GetTotalEnergy() = 0;
    inline static std::string GetBaseDefaultReadName() {return "BondedPotentialBetweenVertices";}
    inline static std::string GetErrorMessage(std::string info) {
        return "---> error: unknown BondedPotentialBetweenVertices type -- \n";
    }
    

};
//---- a class for no box change
class EmptyBonds : public  AbstractBondedPotentialBetweenVertices {
public:
    EmptyBonds(){
        
    }
    ~EmptyBonds(){
        
    }

    inline std::string GetDerivedDefaultReadName() {return "No";}
    inline static std::string GetDefaultReadName() {return "No";}
    void Initialize(){
     
        return;
    }
  /*  double GetBondEnergyOfVertex(vertex *pvertex){
        return 0;
    }*/
    double GetTotalEnergy(){
        return 0;
    }

    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};

#endif
