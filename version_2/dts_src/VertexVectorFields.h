#ifndef VertexVectorFields_H
#define VertexVectorFields_H

#include "SimDef.h"
#include "Vec3D.h"
#include "VectorField.h"

/*
 * @class VertexVectorFields
 * @brief Manages vector fields associated with vertices in a mesh.
 *
 * The VertexVectorFields class is responsible for initializing and managing
 * vector fields that are associated with the vertices of a mesh. This class
 * handles the creation, deletion, and initialization of vector fields based
 * on provided data. The vector fields are represented by instances of the
 * VectorField class.
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */
class MESH;
class VertexVectorFields {

public:
    /*
     * @brief Constructs a VertexVectorFields object with the given ID and inclusion type.
     * @param inctype Pointer to the inclusion type associated with this vector field.
     */
    VertexVectorFields();
    ~VertexVectorFields();

    inline double GetNumberOfVF()               { return m_NoFields; }
    
    
    // this function returns a specific vector field
    inline VectorField * GetVectorField(int vd_layer){
        
        if(vd_layer >= m_NoFields){
            std::cout<<"---> error: this should not be asked \n";
        }
        
    return m_VectorFields[vd_layer];
    }
    
    /*
     * @brief Initializes the vector fields.
     *
     * This function initializes the vector fields with the given number of vector fields
     * and data. It parses the data to create VectorField objects and associates them with
     * the provided mesh.
     *
     * @param no_v The number of vector fields to initialize.
     * @param data A string containing the initialization data for the vector fields.
     * @param pMesh A pointer to the mesh object.
     * @return true if initialization is successful, false otherwise.
     */
    bool Initialize(int no_v, std::string data, MESH *pMesh);
    
    /*
     * @brief Get the vector fields as a formatted string.
     *
     * This function collects all vector fields associated with the vertices
     * and returns them as a formatted string.
     *
     * @return A string representing the vector fields data.
     */
    std::string GetVectorFieldsStream();
    //double CalculateBindingEnergy(vertex *p_vertex);
    double GetBindingEnergy();

    void Copy_VFsBindingEnergy();
    void Reverse_VFsBindingEnergy();

protected:
    std::vector<VectorField *> m_VectorFields;
    int m_NoFields;
    MESH *m_pMesh;

};

#endif // VertexVectorFields_H
