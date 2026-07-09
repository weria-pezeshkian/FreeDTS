#include <stdio.h>
#include "AnalysisCalculations.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "triangle.h"
#include "links.h"



AnalysisCalculations::AnalysisCalculations(State *pState){
    m_pState=pState;
}

void AnalysisCalculations::Calculate(){


    double area_v;
    double area_t=0.0;
    InitializeMemberVariables();

    //Calculating projected area
    if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()){
        m_pProjectedArea=(*(m_pState->GetMesh()->GetBox()))(1)*(*(m_pState->GetMesh()->GetBox()))(0);
    }
    //Calculating total energy
    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()){
        double totalE = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
        m_pState->GetEnergyCalculator()->UpdateTotalEnergy(totalE);
        m_pEnergy=m_pState->GetEnergyCalculator()->GetEnergy();
    }

    //Calculating area and volume
    if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive() || m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
        std::vector<triangle *>& all_tri = m_pState->GetMesh()->GetActiveT();
        for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it) {
             if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()){
                m_pArea += (*it)->GetArea();
             }
             if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()){
                m_pVolume += CalculateSingleTriangleVolume((*it));}
            }
    }



    

    std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();
    //double minHeight=(*all_vertex.begin())->GetVZPos();
    //double maxHeight=(*all_vertex.begin())->GetVZPos();
    double average_height=0.0;
    double inclusion_curvature=0.0;
    double membrane_curvature=0.0;
    //Calculating curvature and thickness
    for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {

        area_v=(*it)->GetArea();
        area_t+=area_v;
        
        //Mean Curvature loop
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature+= area_v*0.5*((*it)->GetP1Curvature() + (*it)->GetP2Curvature());
            m_pMeanCurvature2+= area_v*((*it)->GetP1Curvature() + (*it)->GetP2Curvature())*((*it)->GetP1Curvature() + (*it)->GetP2Curvature())*0.25;

        }
        //Gaussian Curvature loop
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature+= area_v*((*it)->GetP1Curvature()*(*it)->GetP2Curvature());
        }
        //Thickness loop
        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            double z_pos=(*it)->GetVZPos()*area_v;
            average_height+=z_pos; 
        }
    }
    
    if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
        average_height/=area_t;
        for (std::vector<vertex *>::iterator it = all_vertex.begin() ; it != all_vertex.end(); ++it) {
            double z_pos=(*it)->GetVZPos();
            double area_v=(*it)->GetArea();
            m_pThickness+=(z_pos - average_height)*(z_pos - average_height)*area_v;
        }
        m_pThickness=m_pThickness/area_t;
        }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        CalculateInclusion();
        
        if (m_pState->GetAnalysisVariables()->GetNumberInclusionTypes()>1){
            CalculateInclusionNeighbours();
        }
    }




    /*if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature/=area_t;

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature/= area_t;
        }*/

}

void AnalysisCalculations::CalculateInclusionNeighbours(){
    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();

    m_pInclusionNeighbourVector.clear();
    m_pInclusionMeanCurvatureVector.clear();
    m_pInclusionGaussianCurvatureVector.clear();
    m_pInclusionEnergyVector.clear();
    m_pNumberOfInclusionPerTypesVector.clear();
    double area_v;
    int inctype1,inctype2;
    int number_of_neighbours=0;
    int number_of_inclusion_types=m_pState->GetAnalysisVariables()->GetNumberInclusionTypes();
    m_pInclusionNeighbourVector.resize(number_of_inclusion_types, std::vector<double>(number_of_inclusion_types, 0.0));
    m_pInclusionMeanCurvatureVector.resize(number_of_inclusion_types, 0.0);
    m_pInclusionGaussianCurvatureVector.resize(number_of_inclusion_types, 0.0);
    m_pInclusionEnergyVector.resize(number_of_inclusion_types, 0.0);
    m_pNumberOfInclusionPerTypesVector.resize(number_of_inclusion_types, 0);

    for (std::vector<inclusion *>::iterator it = all_inclusions.begin() ; it != all_inclusions.end(); ++it){

        inctype1=(*it)->GetInclusionType()->ITid;
        std::vector<vertex *> neighbour_list= (*it)->Getvertex()->GetVNeighbourVertex();
        for (std::vector<vertex *>::iterator it2 = neighbour_list.begin() ; it2 != neighbour_list.end(); ++it2){
            if ((*it2)->VertexOwnInclusion()){
                inctype2=(*it2)->GetInclusion()->GetInclusionType()->ITid;
                if (inctype1!=inctype2){
                    m_pInclusionNeighbourVector[inctype1-1][inctype2-1]+=1.0;
                    m_pInclusionNeighbourVector[inctype2-1][inctype1-1]+=1.0;
                }
                else{
                    m_pInclusionNeighbourVector[inctype1-1][inctype2-1]+=1.0;
                }
        }
                       
        }

        int inctype1_normalized=inctype1-1;
        m_pInclusionMeanCurvatureVector[inctype1_normalized]+=0.5*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature());
        m_pInclusionGaussianCurvatureVector[inctype1_normalized]+=((*it)->Getvertex()->GetP1Curvature()*(*it)->Getvertex()->GetP2Curvature());
        m_pInclusionEnergyVector[inctype1_normalized]+=(*it)->Getvertex()->GetEnergy();
        m_pNumberOfInclusionPerTypesVector[inctype1_normalized]++;

    }
    
    double divisor=1;
        
        for (size_t i = 0; i < m_pInclusionNeighbourVector.size(); ++i) {
            for (size_t j = 0; j < m_pInclusionNeighbourVector[i].size(); ++j) {
                if (i==j){
                    divisor=static_cast<double>(m_pNumberOfInclusionPerTypesVector[i]);
                }
                else{
                    divisor=static_cast<double>(m_pNumberOfInclusionPerTypesVector[i]+m_pNumberOfInclusionPerTypesVector[j]);
                }
                m_pInclusionNeighbourVector[i][j] =m_pInclusionNeighbourVector[i][j]/divisor; // Divide each element by the divisor
            }
        }
    
    
}

void AnalysisCalculations::CalculateInclusion(){

    std::vector<inclusion *>& all_inclusions = m_pState->GetMesh()->GetInclusion();
    std::vector<vertex *>& all_vertex = m_pState->GetMesh()->GetActiveV();

    double area_v;
    int inctype1,inctype2;

    for (std::vector<inclusion *>::iterator it = all_inclusions.begin() ; it != all_inclusions.end(); ++it){

        //Calculating inclusion neighbours
        std::vector<vertex *> neighbour_list= (*it)->Getvertex()->GetVNeighbourVertex();
        for (std::vector<vertex *>::iterator it2 = neighbour_list.begin() ; it2 != neighbour_list.end(); ++it2){
            if ((*it2)->VertexOwnInclusion()){
                    m_pInclusionNeighbour++;
            }
            
        }


        area_v=(*it)->Getvertex()->GetArea();
        m_pInclusionMeanCurvature+= area_v*0.5*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature());
        m_pInclusionMeanCurvature2+= area_v*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature())*((*it)->Getvertex()->GetP1Curvature() + (*it)->Getvertex()->GetP2Curvature())*0.25;
        m_pInclusionGaussianCurvature+= area_v*((*it)->Getvertex()->GetP1Curvature()*(*it)->Getvertex()->GetP2Curvature());
        m_pInclusionEnergy+=(*it)->Getvertex()->GetEnergy();
        
    }

    m_pInclusionNeighbour/=static_cast<double>(all_inclusions.size());

}


void AnalysisCalculations::InitializeMemberVariables(){

    if (m_pState->GetAnalysisVariables()->GetEnergyCalculationActive()) {
            m_pEnergy=0;
        }
        if (m_pState->GetAnalysisVariables()->GetAreaCalculationActive()) {
            m_pArea=0;
        }
        if (m_pState->GetAnalysisVariables()->GetVolumeCalculationActive()) {
            m_pVolume=0;
        }
        if (m_pState->GetAnalysisVariables()->GetProjectedAreaCalculationActive()) {
            m_pProjectedArea=0;

        }
        if (m_pState->GetAnalysisVariables()->GetMeanCurvatureCalculationActive()) {
            m_pMeanCurvature=0;
            m_pMeanCurvature2=0; // This is used for the mean curvature of inclusions

        }
        if (m_pState->GetAnalysisVariables()->GetGaussianCurvatureCalculationActive()) {
            m_pGaussianCurvature=0;
        }

        if (m_pState->GetAnalysisVariables()->GetThicknessCalculationActive()) {
            m_pThickness=0;
        }

        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
            m_pInclusionNeighbour=0;
            m_pInclusionEnergy=0;
            m_pInclusionMeanCurvature=0;
            m_pInclusionMeanCurvature2=0;
            m_pInclusionGaussianCurvature=0;
        }

}


double AnalysisCalculations::CalculateSingleTriangleVolume(triangle *pTriangle){

    /*if(m_pState->GetMesh()->GetHasCrossedPBC()){
        *(m_pState->GetTimeSeriesLog()) << "---> the system has crossed the PBC while volume is being calculated.";
        *(m_pState->GetTimeSeriesLog()) << " SOLUTION: Restart the simulation and center the system. Also, activate the command for centering the box.";

         exit(0);
    }*/
    
    double T_area = pTriangle->GetArea();
    Vec3D Normal_v = pTriangle->GetNormalVector();
    Vec3D Pos = pTriangle->GetV1()->GetPos();

    // Compute triangle volume
    return T_area * (Vec3D::dot(Normal_v, Pos)) / 3.0;
}

