
#include <complex>
#include <cmath>
#include <unordered_set>
#include <time.h>
#include "RDF.h"
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
RDF::RDF(State* pState, double dr){
    m_pState = pState;
    m_pNumberofInclusions =(m_pState->GetMesh()->GetInclusion()).size();
    Vec3D* box=m_pState->GetMesh()->GetBox();
    m_Lx= (*box)(0);
    m_Ly= (*box)(1);
    double max_distance=std::sqrt(m_Lx * m_Lx + m_Ly * m_Ly);
    double min_distance=0.5;
    m_dr = dr;
    m_pNumberofBins = 0;
    create_distance_vector(min_distance, max_distance, m_dr);
    create_normalization_vector();
    std::cout<<"RDF: maximum distance: " << max_distance << std::endl;
    std::cout<<"RDF: minimum distance: " << min_distance << std::endl;
    std::cout<< "RDF: Number of bins created: " << m_pNumberofBins << std::endl;
}
RDF::~RDF() {
    
}

void RDF::create_distance_vector(double r_min, double r_max, double dr) {
    for (double val = r_min; val <= r_max; val += dr) {
        m_Rvector.push_back(val);
        m_pNumberofBins++;
    }
}
void RDF::create_normalization_vector(){

    for (size_t i = 0; i < m_Rvector.size(); ++i) {
        double r = m_Rvector[i];
        // Calculate the normalization factor for each distance
        double normalizationFactor = 1/ (2 * M_PI * r * m_dr); 
        m_Nvector.push_back(normalizationFactor);
    }
}

void RDF::OpenOutputStreams(bool clearfile) {

    std::string filename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameRDFFile();

    if (!clearfile) {
        m_InclusionCluster.open(filename,std::ios_base::app);  
    }
    else{
        m_InclusionCluster.open(filename);
    }

}

void RDF::CloseOutputStreams() {


    if (m_InclusionCluster.is_open()) {
        m_InclusionCluster.flush(); 
        m_InclusionCluster.close();
     }

}

void RDF::CalculateRDF(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    std::vector<double> histogram(m_pNumberofBins,0);


    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();


    for (inclusion* inc : all_inclusions) {
        double x1 = inc->Getvertex()->GetVXPos();
        double y1 = inc->Getvertex()->GetVYPos();
        for (inclusion* inc2 : all_inclusions) {
            if (inc->GetID() == inc2->GetID()) continue; // Skip self-comparison

            double x2 = inc2->Getvertex()->GetVXPos();
            double y2 = inc2->Getvertex()->GetVYPos();

            // Calculate the distance between the two inclusions
            double dx = x1 - x2;
            double dy = y1 - y2;

            // Apply periodic boundary conditions
            if (dx > m_Lx / 2) dx -= m_Lx;
            if (dx < -m_Lx / 2) dx += m_Lx;
            if (dy > m_Ly / 2) dy -= m_Ly;
            if (dy < -m_Ly / 2) dy += m_Ly;

            double distance = std::sqrt(dx * dx + dy * dy);

            // Find the bin for this distance
            int bin_index = static_cast<int>(distance / m_dr);
            if (bin_index >= 0 && bin_index < m_pNumberofBins) {
                histogram[bin_index] += 1.0; // Increment the count for this bin
            }
        }
        }


    
    // Write the qvector and hqvector to files
    

    if (m_InclusionCluster.is_open()) {
        for (size_t i = 0; i < m_pNumberofBins; ++i) {
            m_InclusionCluster << histogram[i]*m_Nvector[i] << " ";
        }
        m_InclusionCluster << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }


    
}


