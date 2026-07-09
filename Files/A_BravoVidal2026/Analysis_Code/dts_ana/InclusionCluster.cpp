
#include <complex>
#include <cmath>
#include <unordered_set>
#include <time.h>
#include "InclusionCluster.h"
#include "State.h"
#include <algorithm> // For std::sort


using Complex = std::complex<double>;
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Energy of a single vertex
 Energy of a link (when connected vertices has inclusions)
 Energy of the whole system
 */
InclusionCluster::InclusionCluster(State* pState){
    m_pState = pState;
    m_pNumberofInclusions =(m_pState->GetMesh()->GetInclusion()).size();
    

}
InclusionCluster::~InclusionCluster() {
    
}

void InclusionCluster::OpenOutputStreams(bool clearfile) {

    std::string filename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameInclusionClusterFile();

    if (!clearfile) {
        m_InclusionCluster.open(filename,std::ios_base::app);  
    }
    else{
        m_InclusionCluster.open(filename);
    }
}

void InclusionCluster::CloseOutputStreams() {


    if (m_InclusionCluster.is_open()) {
        m_InclusionCluster.flush(); 
        m_InclusionCluster.close();
     }

}

void InclusionCluster::CalculateClusterDistribution(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    std::vector<int> clustervector(m_pNumberofInclusions,0);

    std::unordered_set<inclusion*> visited;

    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();


    for (inclusion* inc : all_inclusions) {
        if (visited.find(inc) != visited.end()) continue;  // Skip already visited inclusions

        // Start DFS
        std::vector<inclusion*> stack = {inc};
        int clusterSize = 0;

        while (!stack.empty()) {
            inclusion* current = stack.back();
            stack.pop_back();

            if (visited.count(current)) continue;
            visited.insert(current);
            clusterSize++;

            // Check neighbors
            std::vector<vertex*> neighbors = current->Getvertex()->GetVNeighbourVertex();
            for (vertex* v : neighbors) {
                if (v->VertexOwnInclusion()) {
                    inclusion* neighborInclusion = v->GetInclusion();
                    if (neighborInclusion && visited.find(neighborInclusion) == visited.end()) {
                        stack.push_back(neighborInclusion);
                    }
                }
            }
        }

        // Store cluster count
        if (clusterSize > 0 && clusterSize < clustervector.size()) {
            clustervector[clusterSize]++;
        }
}


    
    
    // Write the qvector and hqvector to files
    

    if (m_InclusionCluster.is_open()) {

        for (size_t i = 0; i < clustervector.size(); ++i) {
            m_InclusionCluster << clustervector[i] << " ";
        }
        m_InclusionCluster << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }


    
}


