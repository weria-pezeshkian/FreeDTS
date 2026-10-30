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



bool VertexVectorFields::Initialize(int no_v, const std::string& data, MESH* pMesh)
{
    // Initialize the VertexVectorFields with the given number of vector fields and data.


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
    if (no_v <= 0 || pMesh == nullptr) {
        std::cerr << "---> error: invalid mesh pointer, should not happen \n";
        return false;
    }

    m_pMesh = pMesh;

    const auto& all_incType = m_pMesh->GetInclusionType();

    std::vector<std::string> data_str = Nfunction::Split(data);

    constexpr int fields_per_entry = 3;

    if (static_cast<int>(data_str.size()) != fields_per_entry * no_v) {
        std::cerr << "---> error: insufficient or malformed vector field data\n";
        return false;
    }

    // Parse everything first BEFORE modifying object state (strong exception safety)
    struct ParsedField {
        int inc_type_id;
        double x;
        double y;
    };

    std::vector<ParsedField> parsed;
    parsed.reserve(no_v);

    for (int i = 0; i < no_v; ++i) {
        const int base = i * fields_per_entry;

        int inc_type_id = 0;
        double x = 0.0, y = 0.0;

        try {
            inc_type_id = Nfunction::String_to_Int(data_str[base]);
            x = Nfunction::String_to_Double(data_str[base + 1]);
            y = Nfunction::String_to_Double(data_str[base + 2]);
        }
        catch (const std::exception& e) {
            std::cerr << "---> error: failed to parse vector field data: " << e.what() << "\n";
            return false;
        }

        if (inc_type_id < 0 || inc_type_id >= static_cast<int>(all_incType.size())) {
            std::cerr << "\033[1;31m [ERROR]-->\033[0m : vector field type does not exist \n";
            return false;
        }

        parsed.push_back({inc_type_id, x, y});
    }

    // Now safely replace internal state
    std::vector<std::unique_ptr<VectorField>> newFields;
    newFields.reserve(no_v);

    for (int i = 0; i < no_v; ++i) {
        const auto& f = parsed[i];
        newFields.push_back(std::make_unique<VectorField>(
            i,
            all_incType[f.inc_type_id],
            f.x,
            f.y
        ));
    }

    // Commit only after success
    m_VectorFields.clear();
    for (auto& f : newFields) {
        m_VectorFields.push_back(f.release());
    }

    m_NoFields = no_v;

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
