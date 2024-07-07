#include <stdio.h>
#include "VertexVectorFields.h"
#include "MESH.h"

// Constructor: Initialize the number of fields to zero and clear the vector fields.
VertexVectorFields::VertexVectorFields() {
    m_NoFields = 0;
    m_VectorFields.clear();
}

// Destructor: Clean up dynamically allocated VectorField objects.
VertexVectorFields::~VertexVectorFields() {
    for (size_t i = 0; i < m_VectorFields.size(); ++i) {
        delete m_VectorFields[i];
    }
    m_VectorFields.clear(); // Clear the vector to remove dangling pointers.
}

// Initialize the VertexVectorFields with the given number of vector fields and data.
bool VertexVectorFields::Initialize(int no_v, std::string data, MESH *pMesh) {
    /**
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
    m_pMesh = pMesh;
    std::vector<InclusionType*> all_incType = m_pMesh->GetInclusionType();

    std::vector<std::string> data_str = Nfunction::Split(data);
    m_NoFields = no_v;

    // Check if the data provided matches the expected format.
    if (data_str.size() != 3 * m_NoFields) {
        std::cout << "---> error, info on vector field is not enough \n";
        return false;
    }

    // Clear existing vector fields if any.
    for (size_t i = 0; i < m_VectorFields.size(); ++i) {
        delete m_VectorFields[i];
    }
    m_VectorFields.clear();
    m_VectorFields.reserve(m_NoFields); // Reserve space for new vector fields.

    // Create new VectorField objects and add them to m_VectorFields.
    for (int i = 0; i < no_v; ++i) {
        int inc_type_id = Nfunction::String_to_Int(data_str[i * 3]);

        if (inc_type_id >= static_cast<int>(all_incType.size())) {
            std::cout << "---> vector field: inc type id too high \n";
            return false;
        }

        double x = Nfunction::String_to_Double(data_str[i * 3 + 1]);
        double y = Nfunction::String_to_Double(data_str[i * 3 + 2]);
        InclusionType* inc_type = all_incType[inc_type_id];

        // Allocate a new VectorField object and add it to the vector.
        m_VectorFields.push_back(new VectorField(i, inc_type, x, y));
    }

    return true;
}
std::string VertexVectorFields::GetVectorFieldsStream(){
    /*
     * @brief Get the vector fields as a formatted string.
     *
     * This function collects all vector fields associated with the vertices
     * and returns them as a formatted string.
     *
     * @return A string representing the vector fields data.
     */
    std::string str_data;
    str_data.reserve(m_VectorFields.size() * 30); // Reserve memory to reduce reallocations

    for (std::vector<VectorField*>::iterator it = m_VectorFields.begin(); it != m_VectorFields.end(); ++it) {
        double x = (*it)->GetLDirection()(0);
        double y = (*it)->GetLDirection()(1);
        int tid = (*it)->GetInclusionType()->ITid;

        str_data += "   " + Nfunction::D2S(tid) + "   " + Nfunction::D2S(x) + "   " + Nfunction::D2S(y);
    }
    
    return str_data;
}
/*double VertexVectorFields::CalculateBindingEnergy(vertex *p_vertex){
    double T_en = 0;
    
    for (std::vector<VectorField*>::iterator it = m_VectorFields.begin(); it != m_VectorFields.end(); ++it) {
        double en = (*it)->CalculateMembraneBindingEnergy(p_vertex);
        (*it)->UpdateMembraneBindingEnergy(en);
        T_en += en;
    }
    
    return T_en;
}*/
double VertexVectorFields::GetBindingEnergy(){
 
    if(m_NoFields == 0 ){
        return 0;
    }
        
    double en = 0;
    for (std::vector<VectorField*>::iterator it = m_VectorFields.begin(); it != m_VectorFields.end(); ++it) {
        en += (*it)->GetMembraneBindingEnergy();
    }
    
    return en;
}
void VertexVectorFields::Copy_VFsBindingEnergy(){
    
    if(m_NoFields == 0 ){
        return;
    }
    for (std::vector<VectorField*>::iterator it = m_VectorFields.begin(); it != m_VectorFields.end(); ++it) {
        (*it)->Copy_BindingEnergy();
    }
    return;
}
void VertexVectorFields::Reverse_VFsBindingEnergy(){
    
    if(m_NoFields == 0 ){
        return;
    }
    for (std::vector<VectorField*>::iterator it = m_VectorFields.begin(); it != m_VectorFields.end(); ++it) {
        (*it)->Reverse_BindingEnergy();
    }
    return;
}
