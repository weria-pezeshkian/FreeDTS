#include "ConstantExternalFieldOnVectorFields.h"

// E = -k*(Field*Inc_direction)^2 --> we could make it k1*(Field*Inc_direction)-k2*(Field*Inc_direction)^2

ConstantExternalFieldOnVectorFields::ConstantExternalFieldOnVectorFields(std::string data_stream) {

    std::vector<std::string> data = Nfunction::Split(data_stream);
    if(data.size() == 0 || data.size()%4 != 0 ){
        std::cout<<"---> error: data provided for ConstantExternalFieldOnVectorFields is not enough \n";
        std::cout<<"-------- the data should be 4x doubles \n";

    }
    for (int i = 0; i<data.size(); i=i+4){
        double f = Nfunction::String_to_Double(data[i]);
        double x = Nfunction::String_to_Double(data[i+1]);
        double y = Nfunction::String_to_Double(data[i+2]);
        double z = Nfunction::String_to_Double(data[i+3]);

        m_vFieldStrength.push_back(f);
        double norm = sqrt(x*x + y*y + z*z);
        Vec3D tem(x, y, z);
        tem.normalize(); // Scale the vector by the reciprocal of its magnitude to normalize it
        m_vFieldDirection.push_back(tem);
    }
}
ConstantExternalFieldOnVectorFields::~ConstantExternalFieldOnVectorFields() {
    
}
double ConstantExternalFieldOnVectorFields::GetCouplingEnergy(int layer, VectorField* p_vf, vertex *pvertex) {
    // Check if the vertex owns an inclusion

    if(layer>= m_vFieldStrength.size()){
        std::cout<<"--->  error 77373 \n";
    }
    
    if(m_vFieldStrength[layer] == 0){
        return 0; // return zero coupling energy
    }

    // Get the local direction of the inclusion
    Vec3D LD = p_vf->GetLDirection();

    // Transform the local direction to global direction using the transfer matrix
    Vec3D GD = pvertex->GetL2GTransferMatrix() * LD;

    // Calculate the cosine of the angle between the global direction and the field direction
    double Cangle = GD.dot(GD,m_vFieldDirection[layer]);

    // Calculate and return the coupling energy
    // Note: The negative sign indicates that the field tends to minimize the angle with the field direction
    return -m_vFieldStrength[layer] * Cangle * Cangle;
}
std::string ConstantExternalFieldOnVectorFields::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    for (int i = 0; i<m_vFieldStrength.size(); i++) {
        
        state += "  "+ Nfunction::D2S(m_vFieldStrength[i]);
        state += "  "+ Nfunction::D2S(m_vFieldDirection[i](0));
        state += "  "+ Nfunction::D2S(m_vFieldDirection[i](1));
        state += "  "+ Nfunction::D2S(m_vFieldDirection[i](2));

    }
    return state;
}
