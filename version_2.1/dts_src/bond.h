#ifndef BOND_H
#define BOND_H

#include "SimDef.h"
#include "triangle.h"
#include "Vec3D.h"
#include "Tensor2.h"

/*
 * -----------------------------------------------------------------------------
 * Bond Class - Simulation Component
 * -----------------------------------------------------------------------------
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Date: December 2024
 *
 * Description:
 * The `bond` class represents a bond in a  simulation. It models
 * the connection between two vertices with properties such as stiffness
 * (k) and rest length (l0). The class provides methods to calculate
 * bond energy and update bond parameters dynamically during the simulation.
 *  E = k/2(l-l0)^2
 *
 * -----------------------------------------------------------------------------
 * License:
 * Copyright (c) Weria Pezeshkian, 2024. All rights reserved.
 * -----------------------------------------------------------------------------
 */

class vertex;

class bond {
public:
    // Constructors and Destructor
    bond(int id, vertex* v1, vertex* v2, double k, double l0);
    explicit bond(int id);
    ~bond();

    // Accessor Methods
    int GetID() const { return m_ID; }
    vertex* GetV1() const { return m_V1; }
    vertex* GetV2() const { return m_V2; }

    // Public Member Functions
    void UpdateV(vertex* v1, vertex* v2);
    void UpdateBondK(double k);
    void UpdateBondL(double l0);
    double CalculateEnergy() const;

private:
    // Member Variables
    int m_ID;
    vertex* m_V1;
    vertex* m_V2;
    double m_K;
    double m_L0;
    Vec3D *m_pBox; 

};

#endif // BOND_H
