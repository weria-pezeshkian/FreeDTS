#include <cmath>
#include <time.h>
#include "InclusionTracking.h"
#include "State.h"
#include <algorithm> // For std::sort


InclusionTracking::InclusionTracking(State* pState){
    m_pState = pState;
    OpenOutputStreams(true);

}
InclusionTracking::~InclusionTracking() {
    
}

void InclusionTracking::OpenOutputStreams(bool clearfile) {

    //Change this line to include the type of file that you want
    std::string hqfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameInclusionTrackingFile();


    if (!clearfile) {
        m_HQVector.open(hqfilename,std::ios_base::app);
       
    }
    else{
        m_HQVector.open(hqfilename);
    }
}

void InclusionTracking::CloseOutputStreams() {


    if (m_HQVector.is_open()) {
        m_HQVector.flush(); 
        m_HQVector.close();
     }

}





void InclusionTracking::OutputInclusionInfo(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    //Just write out everything related to each inclusion
    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();
    double area_v, m_pInclusionMeanCurvature=0.0, m_pInclusionEnergy=0.0;
    int inctype1=0;
    int idx=0;


    if (m_HQVector.is_open()){
    for (std::vector<inclusion *>::iterator it = all_inclusions.begin() ; it != all_inclusions.end(); ++it){
        idx=(*it)->GetID();
        area_v=(*it)->Getvertex()->GetArea();
        m_pInclusionMeanCurvature=0.5*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature());
        m_pInclusionEnergy=(*it)->Getvertex()->GetEnergy();
        inctype1=(*it)->GetInclusionType()->ITid;
        m_HQVector<<inctype1<<","<<area_v<<","<<m_pInclusionMeanCurvature<<","<<m_pInclusionEnergy<<" ";

    }
    m_HQVector<<"\n";
    } else {
        std::cerr << "Unable to open inclusiontracking.txt for writing." << std::endl;
    }



    //if (m_HQVector.is_open()) {
    //    for (size_t i = 0; i < hqvector.size(); ++i) {
    //        m_HQVector << hqvector[i] << " ";
    //    }
    //    m_HQVector << "\n";
    //} else {
    //    std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    //}

    
    
}


