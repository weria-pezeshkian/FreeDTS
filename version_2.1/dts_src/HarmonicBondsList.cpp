#include "HarmonicBondsList.h"
#include "State.h"

// E = -k*(Field*Inc_direction)^2 --> we could make it k1*(Field*Inc_direction)-k2*(Field*Inc_direction)^2

HarmonicBondsList::HarmonicBondsList(State *pState, std::string filename) {
    m_FileName = filename;
    m_pState = pState;
}
HarmonicBondsList::~HarmonicBondsList()
{
    
}
void HarmonicBondsList::Initialize() {
    // Open the bond file
    std::ifstream input(m_FileName);
    if (!input) {
        std::cerr << "Error: Failed to open the input file: " << m_FileName << std::endl;
        exit(EXIT_FAILURE);
    }

    // Get active vertices
    const std::vector<vertex*>& activeVertices = m_pState->GetMesh()->GetActiveV();
    const size_t numVertices = activeVertices.size();

    std::string line;
    int bondId = 0;

    // Read the file line by line
    while (std::getline(input, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == ';') {
            continue;
        }

        // Split the line into components
        std::vector<std::string> tokens = Nfunction::Split(line);
        if (tokens.size() < 3) {
            std::cerr << "Error: Malformed line in bond list: " << line << std::endl;
            exit(EXIT_FAILURE);
        }

        // Parse bond information
        int id1 = Nfunction::String_to_Int(tokens[0]);
        int id2 = Nfunction::String_to_Int(tokens[1]);
        double k = Nfunction::String_to_Double(tokens[2]);
        double l0 = (tokens.size() > 3) ? Nfunction::String_to_Double(tokens[3]) : 0.0; // Default l0 if missing

        // Validate vertex indices
        if (id1 >= numVertices || id2 >= numVertices) {
            std::cerr << "Error: Vertex index out of range in bond list: " << line << std::endl;
            exit(EXIT_FAILURE);
        }

        // Create and store the bond
        bond b(bondId++, activeVertices[id1], activeVertices[id2], k, l0);
        m_AllBonds.push_back(b);
    }

    // Associate bonds with their vertices
    for (bond& b : m_AllBonds) {
        b.GetV1()->AddBondToList(&b);
        b.GetV2()->AddBondToList(&b);
    }

}
std::string HarmonicBondsList::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+m_FileName;
    return state;
}
double HarmonicBondsList::GetTotalEnergy(){
    
    double en = 0;
    
    for (std::vector<bond>::iterator it = m_AllBonds.begin() ; it != m_AllBonds.end(); ++it){
        en += it->CalculateEnergy();
    }
    
    return en;
}
