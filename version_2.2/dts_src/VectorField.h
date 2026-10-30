#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include "SimDef.h"
#include "Vec3D.h"
#include "InclusionType.h"

/*
 * @brief The VectorField class represents a vector field on a vertex.
 * A vector field is similar to an inclusion but does not jump and each vertex could have multiple 
 * This class models a vector field with local and global directional properties.
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */
class vertex;
class VectorField {

public:
    /**
     * @brief Constructs a VectorField object with the given ID and inclusion type.
     *
     * @param inctype Pointer to the inclusion type associated with this vector field.
     */
    VectorField();
    VectorField(int layer, InclusionType* inctype, double x, double y);
    ~VectorField();

    inline Vec3D GetLDirection() { return m_LDirection; } //Gets the direction of the vector field in the local frame.
    inline Vec3D GetGDirection() { return m_GDirection; } //Gets the direction of the vector field in the local frame.

    inline InclusionType* GetInclusionType() { return m_IncType; } //Gets the direction of the vector field in the local frame.
    inline double GetMembraneBindingEnergy() { return m_MembraneBindingEnergy; } 
    inline int GetLayer() { return m_Layer; }

    void UpdateLocalDirection(const Vec3D& lo_dir); //Updates the direction of the vector field in the local frame.
    void UpdateGlobalDirection(const Vec3D& lg_dir); //Updates the direction of the vector field in the global frame.
    //double CalculateMembraneBindingEnergy(vertex *p_vertex);
    void UpdateMembraneBindingEnergy(const double &en);
    void Add2MembraneBindingEnergy(const double &en);
    
    void Copy_BindingEnergy(){
        m_OldMembraneBindingEnergy = m_MembraneBindingEnergy;
    };
    void Reverse_BindingEnergy(){
        m_MembraneBindingEnergy = m_OldMembraneBindingEnergy;
    };
    void Copy_Direction(){
        m_OldLDirection = m_LDirection;
    };
    void Reverse_Direction(){
        m_LDirection = m_OldLDirection;
    };
    
private:
    Vec3D m_LDirection;         ///< Direction of the vector field in the local frame.
    Vec3D m_OldLDirection;         // a copy for local direction
    Vec3D m_GDirection;         ///< Direction of the vector field in the global frame.
    InclusionType* m_IncType;   ///< Pointer to the inclusion type associated with this vector field.
    double m_MembraneBindingEnergy;
    int m_Layer;
    double m_OldMembraneBindingEnergy;

};

#endif // VECTORFIELD_H
