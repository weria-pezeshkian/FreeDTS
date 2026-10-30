
#include <stdio.h>
#include "VectorField.h"
#include "vertex.h"

VectorField::VectorField()  {

}
VectorField::VectorField(int layer, InclusionType *inctype, double x, double y) : m_IncType(inctype), m_Layer(layer) {
    
    m_LDirection(0) = x;
    m_LDirection(1) = y;
    m_LDirection(2) = 0;
    m_LDirection.normalize();
    m_MembraneBindingEnergy = 0;
}
VectorField::~VectorField() {
    
}
void VectorField::UpdateLocalDirection(const Vec3D & lo_dir) {
    m_LDirection = lo_dir;
    return;
}
void VectorField::UpdateGlobalDirection(const Vec3D & lg_dir) {
    m_GDirection = lg_dir;
    return;
}
void VectorField::UpdateMembraneBindingEnergy(const double &en){
    
    m_MembraneBindingEnergy = en;
    return;
}
void VectorField::Add2MembraneBindingEnergy(const double &en){
    
    m_MembraneBindingEnergy += en;
    return;
}


