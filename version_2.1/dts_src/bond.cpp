#include "bond.h"
#include "vertex.h"

// Constructor with initialization of bond parameters
bond::bond(int id, vertex* v1, vertex* v2, double k, double l0)
    : m_ID(id), m_V1(v1), m_V2(v2), m_K(k / 2), m_L0(l0) {
        
       m_pBox = m_V1->GetBox();
    }

// Constructor for uninitialized bond
bond::bond(int id)
    : m_ID(id), m_V1(nullptr), m_V2(nullptr), m_pBox(nullptr), m_K(0.0), m_L0(0.0) {}

// Destructor
bond::~bond() = default;

// Update vertices associated with the bond
void bond::UpdateV(vertex* v1, vertex* v2) {
    m_V1 = v1;
    m_V2 = v2;
}

// Update the bond stiffness
void bond::UpdateBondK(double k) {
    m_K = k / 2;
}

// Update the bond's rest length
void bond::UpdateBondL(double l0) {
    m_L0 = l0;
}

// Calculate the bond's energy
double bond::CalculateEnergy() const {
    Vec3D displacement = m_V1->GetPos() - m_V2->GetPos();

    for (int i=0;i<3;i++){
        if(fabs(displacement(i))>(*m_pBox)(i)/2.0)
        {
            if(displacement(i)<0)
                displacement(i) = (*m_pBox)(i) + displacement(i);
            else if(displacement(i)>0)
                displacement(i) =  displacement(i) - (*m_pBox)(i);
        }
    }

    double currentLength = displacement.norm();
    double energy = m_K * (currentLength - m_L0) * (currentLength - m_L0);
    return energy;
}
