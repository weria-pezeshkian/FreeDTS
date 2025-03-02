#include <stdio.h>
#include "VertexHarmonicBounds.h"
#include "MESH.h"

VertexHarmonicBounds::VertexHarmonicBounds() {
    m_pVertexBond.clear();
}

VertexHarmonicBounds::~VertexHarmonicBounds() {

}
void VertexHarmonicBounds::AddBondToList(bond *b){
    
    m_pVertexBond.push_back(b);
    return;
}
double VertexHarmonicBounds::GetBondEnergyOfVertex(){
    
    double en = 0;
    
    for (std::vector<bond*>::iterator it = m_pVertexBond.begin() ; it != m_pVertexBond.end(); ++it){
        en += (*it)->CalculateEnergy();
    }
    
    return en;
}
