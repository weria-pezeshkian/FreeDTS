#include <cmath>
#include <time.h>
#include "FluctuationSpectrum.h"
#include "State.h"
#include <algorithm> // For std::sort



struct FourierResult {
    std::vector<Complex> sum_p;
    std::vector<Complex> sum_h;
};
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Energy of a single vertex
 Energy of a link (when connected vertices has inclusions)
 Energy of the whole system
 */
FluctuationSpectrum::FluctuationSpectrum(State* pState,int nx,int ny,int type){
    m_pState = pState;
    m_Nx = nx;
    m_Ny = ny;


    m_type=type;
    if(m_type==0)
        {GenerateVectorOrder();
            m_SpectrumSize=m_VectorOrder.size();
            GenerateZeroAndNonZeroVectorOrder();}
    else if(m_type==1 or m_type==2)
        {GenerateVectorOrder();
        m_SpectrumSize=m_VectorOrder.size();
            GenerateZeroAndNonZeroVectorOrderIndividual();}
    else if(m_type==3){
        m_lmax=nx;
        m_SpectrumSize=(m_lmax+1)*(m_lmax+1);
        GenerateVectorOrderSpherical();
        std::cout<<m_SpectrumSize<<std::endl;
        std::cout<<"Vector Order size: "<<m_VectorOrder.size()<<std::endl;
        }

    

    


}
FluctuationSpectrum::~FluctuationSpectrum() {
    
}

void FluctuationSpectrum::OpenOutputStreams(bool clearfile) {

    std::string hqfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameHQVectorFile();
    std::string qfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameQVectorFile();
    
    std::string hpqfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameHPQVectorFile();
    std::string pqfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNamePQVectorFile();

    if (!clearfile) {
        m_QVector.open(qfilename,std::ios_base::app);
        m_HQVector.open(hqfilename,std::ios_base::app);

        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
            m_HPQVector.open(hpqfilename,std::ios_base::app);
            m_PQVector.open(pqfilename,std::ios_base::app);
        }
        
        
    }
    else{
        m_QVector.open(qfilename);
        m_HQVector.open(hqfilename);
        if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
            m_HPQVector.open(hpqfilename);
            m_PQVector.open(pqfilename);
        }
    }
}

void FluctuationSpectrum::CloseOutputStreams() {


    if (m_QVector.is_open()) {
        m_QVector.flush(); 
        m_QVector.close();
     }

    if (m_HQVector.is_open()) {
        m_HQVector.flush(); 
        m_HQVector.close();
     }

     if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
        if (m_PQVector.is_open()) {
            m_PQVector.flush(); 
            m_PQVector.close();
         }
         if (m_HPQVector.is_open()) {
            m_HPQVector.flush(); 
            m_HPQVector.close();
         }
    }

}

void FluctuationSpectrum::GenerateVectorOrder(){
    std::vector<int> input(2,0);
    for(int i = 0; i < m_Ny; i++){
        for(int j = 0; j <= i; j++){
            input[1]=j;
            input[0]=i;
            m_VectorOrder.push_back(input);
        }
    }

    auto distance = [](const std::vector<int>& v) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1]);
    };

    // Sort m_VectorOrder based on the calculated distance
    std::sort(m_VectorOrder.begin(), m_VectorOrder.end(),
              [&distance](const std::vector<int>& a, const std::vector<int>& b) {
                  return distance(a) < distance(b);
              });
}

void FluctuationSpectrum::GenerateVectorOrderSpherical(){
    std::vector<int> input(2,0);
    int count=0;
    for(int l = 0; l <= m_lmax; l++){
        for(int m = -l; m <= l; m++){
            input[0] = l;  // l value
            input[1] = m;  // m value
            m_VectorOrder.push_back(input);
            count++;
        }
    }
    std::cout<<"Generated "<<count<<" spherical harmonics up to lmax="<<m_lmax<<std::endl;
}

void FluctuationSpectrum::GenerateZeroAndNonZeroVectorOrder(){

    for (size_t i = 0; i < m_VectorOrder.size(); ++i) {
        std::vector<std::vector<int>> MatrixOrder(4, std::vector<int>(2, 0));


        if(m_VectorOrder[i][0]==m_VectorOrder[i][1]){
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][0], -m_VectorOrder[i][1]};
            m_MatrixOrder.push_back(MatrixOrder);}
        else if (m_VectorOrder[i][0]==0 || m_VectorOrder[i][1]==0){
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][1], m_VectorOrder[i][0]};
            m_MatrixOrder.push_back(MatrixOrder);
        }
        
        else{
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][1], m_VectorOrder[i][0]};
            MatrixOrder[2] = {m_VectorOrder[i][0], -m_VectorOrder[i][1]};
            MatrixOrder[3] = {m_VectorOrder[i][1], -m_VectorOrder[i][0]};
            m_MatrixOrder.push_back(MatrixOrder);
        }

    }
}

void FluctuationSpectrum::GenerateZeroAndNonZeroVectorOrderIndividual(){
    m_MatrixOrderIndividual.push_back({0,0});
    for (size_t i = 1; i < m_VectorOrder.size(); ++i) {
        if(m_VectorOrder[i][0]==m_VectorOrder[i][1]){
            m_MatrixOrderIndividual.push_back(m_VectorOrder[i]);
            m_MatrixOrderIndividual.push_back({m_VectorOrder[i][0], -m_VectorOrder[i][1]});
        }
        else if (m_VectorOrder[i][0]==0 || m_VectorOrder[i][1]==0){
            m_MatrixOrderIndividual.push_back(m_VectorOrder[i]);
            m_MatrixOrderIndividual.push_back({m_VectorOrder[i][1], m_VectorOrder[i][0]});
        }
        
        else{
            m_MatrixOrderIndividual.push_back(m_VectorOrder[i]);
            m_MatrixOrderIndividual.push_back({m_VectorOrder[i][1], m_VectorOrder[i][0]});
            m_MatrixOrderIndividual.push_back({m_VectorOrder[i][0], -m_VectorOrder[i][1]});
            m_MatrixOrderIndividual.push_back({m_VectorOrder[i][1], -m_VectorOrder[i][0]});
        }

    }
}

std::vector<double> FluctuationSpectrum::FourierTransformInclusion(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_h(qvector.size(),Complex(0,0));
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    bool r;
    double result_h=0;
    double result_hp=0;
    double result_p=0;
    std::vector<double> final_result(3,0);


    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            r = (*it)->VertexOwnInclusion() ? 1.0 : 0.0;
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(r,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            sum_h[i]=sum_h[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    for (size_t i = 0; i < qvector.size(); ++i) {
        result_p+=(sum_p[i]*std::conj(sum_p[i])).real();
        result_h+=(sum_h[i]*std::conj(sum_h[i])).real();
        result_hp+=(sum_h[i]*std::conj(sum_p[i])).real();
        //std::cout<<"Fourier result: "<<result<<std::endl;
    }

    result_p=result_p/static_cast<double>(qvector.size());
    result_h=result_h/static_cast<double>(qvector.size());
    result_hp=result_hp/static_cast<double>(qvector.size());

    final_result[0]=result_h;
    final_result[1]=result_hp;
    final_result[2]=result_p;

    return final_result;
}

std::vector<std::vector<double>> FluctuationSpectrum::FourierTransformInclusionIndividual(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_h(qvector.size(),Complex(0,0));
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    bool r;
    double result_h=0;
    double result_hp=0;
    double result_p=0;
    std::vector<std::vector<double>> final_result(3, std::vector<double>(qvector.size(),0.0));


    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            r = (*it)->VertexOwnInclusion() ? 1.0 : 0.0;
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(r,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            sum_h[i]=sum_h[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    for (size_t i = 0; i < qvector.size(); ++i) {
        final_result[0][i]=(sum_p[i]*std::conj(sum_p[i])).real();
        final_result[1][i]=(sum_h[i]*std::conj(sum_h[i])).real();
        final_result[2][i]=(sum_h[i]*std::conj(sum_p[i])).real();
        //std::cout<<"Fourier result: "<<result<<std::endl;
    }



    return final_result;
}

FourierResult FluctuationSpectrum::FourierTransformInclusionIndividualComplex(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    FourierResult result;
    result.sum_h = std::vector<Complex>(qvector.size(),Complex(0,0));
    result.sum_p = std::vector<Complex>(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    bool r;
    double result_h=0;
    double result_hp=0;
    double result_p=0;
    int upper=0;
    int lower=0;



    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();

            

            
            r = (*it)->VertexOwnInclusion() ? 1.0 : 0.0;
            for (size_t i = 0; i < qvector.size(); ++i) {
            result.sum_p[i]=result.sum_p[i]+Complex(r,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            result.sum_h[i]=result.sum_h[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    return result;
}

double FluctuationSpectrum::FourierTransformNoInclusion(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    double result=0;

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    for (size_t i = 0; i < qvector.size(); ++i) {
        result+=(sum_p[i]*std::conj(sum_p[i])).real();
        //std::cout<<"Fourier result: "<<result<<std::endl;
    }

    result=result/static_cast<double>(qvector.size());

    return result;
}

std::vector<double> FluctuationSpectrum::FourierTransformNoInclusionIndividual(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    std::vector<double> result(qvector.size(), 0.0);

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    for (size_t i = 0; i < qvector.size(); ++i) {
        result[i]=(sum_p[i]*std::conj(sum_p[i])).real();
        //std::cout<<"Fourier result: "<<result<<std::endl;
    }


    return result;
}

std::vector<Complex> FluctuationSpectrum::FourierTransformNoInclusionIndividualComplex(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }



    return sum_p;
}

FourierResult FluctuationSpectrum::FourierTransformInclusionSphericalComplex(std::vector<std::vector<int>> qvector){
        const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
        FourierResult result;
        result.sum_h = std::vector<Complex>(m_SpectrumSize,Complex(0,0));
        result.sum_p = std::vector<Complex>(m_SpectrumSize,Complex(0,0));

        double x;
        double y;
        double z;
        bool r;
        double result_h=0;
        double result_hp=0;
        double result_p=0;

        double ui;
        double phi;
        double theta;
        double ri;

        for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
                x=(*it)->GetVXPos()-m_CenterOfMass[0];
                y=(*it)->GetVYPos()-m_CenterOfMass[1];
                z=(*it)->GetVZPos()-m_CenterOfMass[2];
                r = (*it)->VertexOwnInclusion() ? 1.0 : 0.0;

                ri=std::sqrt(x*x+y*y+z*z);
                ui=(ri-m_radius0)/m_radius0;
                phi=std::atan2(y,x);
                theta=std::acos(z/ri);

                
                for (size_t i = 0; i < qvector.size(); ++i) {
                    double ylm_real, ylm_imag;
                    SphericalHarmonic(qvector[i][0], qvector[i][1], theta, phi, ylm_real, ylm_imag);
                    result.sum_h[i]=result.sum_h[i]+Complex(ui,0)*Complex(ylm_real, ylm_imag);
                    result.sum_p[i]=result.sum_p[i]+Complex(r,0)*Complex(ylm_real, ylm_imag);
                }

        }

    return result;
}

void FluctuationSpectrum::AverageHeight(){
    
    double average_height=0;
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    double number_vertex=0;
    double Lz=(*(m_pState->GetMesh()->GetBox()))(2);
    double z;
    int upper=0;
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            z=(*it)->GetVZPos();
            average_height+=z;
            number_vertex+=1;
    }

    average_height/=number_vertex;
    m_AverageHeight=average_height;
}

void FluctuationSpectrum::CenterOfMass(){
    double center_x=0;
    double center_y=0;
    double center_z=0;
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    double number_vertex=0;

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            double x=(*it)->GetVXPos();
            double y=(*it)->GetVYPos();
            double z=(*it)->GetVZPos();
            center_x+=x;
            center_y+=y;
            center_z+=z;
            number_vertex+=1;
    }

    center_x/=number_vertex;
    center_y/=number_vertex;
    center_z/=number_vertex;

    m_CenterOfMass={center_x,center_y,center_z};
}

void FluctuationSpectrum::CalculateSpectrum(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    double Lx=(*(m_pState->GetMesh()->GetBox()))(0);
    double Ly=(*(m_pState->GetMesh()->GetBox()))(1);
    AverageHeight();

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        const std::vector<vertex*>& pAllVertices = m_pState->GetMesh()->GetActiveV();
        int number_of_vertices = pAllVertices.size();
        int number_of_inclusions = (m_pState->GetMesh()->GetInclusion()).size(); // Assuming GetInclusions() returns the binclusion vector
        m_AverageInclusionDensity = static_cast<double>(number_of_inclusions) / number_of_vertices;
    }

    std::vector<double> qvector(m_VectorOrder.size(),0);
    std::vector<double> hqvector(m_VectorOrder.size(),0);
    std::vector<double> hpqvector(m_VectorOrder.size(),0);
    std::vector<double> pqvector(m_VectorOrder.size(),0);
    std::vector<double> container(3,0);

    for (size_t i = 0; i < m_MatrixOrder.size(); ++i) {

        std::vector<double> q={m_VectorOrder[i][0]*2*PI/Lx,m_VectorOrder[i][1]*2*PI/Ly};
        qvector[i]=std::sqrt(std::pow(q[0], 2) + std::pow(q[1], 2));

        if((m_VectorOrder[i][0]==m_VectorOrder[i][1]) || (m_VectorOrder[i][0]==0 || m_VectorOrder[i][1]==0)){
            std::vector<std::vector<double>> qmatrix(2,std::vector<double>(2,0));
            for (size_t j = 0; j < 2; ++j) {
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }
            if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
                container=FourierTransformInclusion(qmatrix);
                hqvector[i]=container[0];
                pqvector[i]=container[2];
                hpqvector[i]=container[1];}
            else{
                hqvector[i]=FourierTransformNoInclusion(qmatrix);
            }
        }   
        else{
            std::vector<std::vector<double>> qmatrix(4,std::vector<double>(2,0));
            for (size_t j = 0; j < 4; ++j) {
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }

            if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
                container=FourierTransformInclusion(qmatrix);
                hqvector[i]=container[0];
                pqvector[i]=container[2];
                hpqvector[i]=container[1];}
            else{
                hqvector[i]=FourierTransformNoInclusion(qmatrix);
            } 
            }

    }

    
    // Write the qvector and hqvector to files
    

    if (m_QVector.is_open()) {
        for (size_t i = 0; i < qvector.size(); ++i) {
            m_QVector << qvector[i] << " ";
        }
        m_QVector << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }

    if (m_HQVector.is_open()) {
        for (size_t i = 0; i < hqvector.size(); ++i) {
            m_HQVector << hqvector[i] << " ";
        }
        m_HQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
    if (m_HPQVector.is_open()) {
        for (size_t i = 0; i < hpqvector.size(); ++i) {
            m_HPQVector << hpqvector[i] << " ";
        }
        m_HPQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }

    if (m_PQVector.is_open()) {
        for (size_t i = 0; i < pqvector.size(); ++i) {
            m_PQVector << pqvector[i] << " ";
        }
        m_PQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }   
}
}

void FluctuationSpectrum::CalculateSpectrumIndividual(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    double Lx=(*(m_pState->GetMesh()->GetBox()))(0);
    double Ly=(*(m_pState->GetMesh()->GetBox()))(1);
    AverageHeight();

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        const std::vector<vertex*>& pAllVertices = m_pState->GetMesh()->GetActiveV();
        int number_of_vertices = pAllVertices.size();
        int number_of_inclusions = (m_pState->GetMesh()->GetInclusion()).size(); // Assuming GetInclusions() returns the binclusion vector
        m_AverageInclusionDensity = static_cast<double>(number_of_inclusions) / number_of_vertices;
    }

    std::vector<std::vector<double>> qvector(m_MatrixOrderIndividual.size(), std::vector<double>(2, 0));
    std::vector<double> hqvector(m_MatrixOrderIndividual.size(),0);
    std::vector<double> hpqvector(m_MatrixOrderIndividual.size(),0);
    std::vector<double> pqvector(m_MatrixOrderIndividual.size(),0);

    for (size_t i = 0; i < m_MatrixOrderIndividual.size(); ++i) {
        std::vector<double> q={m_MatrixOrderIndividual[i][0]*2*PI/Lx,m_MatrixOrderIndividual[i][1]*2*PI/Ly};
        qvector[i]=q;
    }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
        std::vector<std::vector<double>> final_result=FourierTransformInclusionIndividual(qvector);
        hqvector = final_result[1];
        hpqvector = final_result[2];
        pqvector = final_result[0];
    }
    else{
        std::vector<double> final_result=FourierTransformNoInclusionIndividual(qvector);
        hqvector = final_result;
    }


    

    
    // Write the qvector and hqvector to files
    

    if (m_QVector.is_open()) {
        for (size_t i = 0; i < qvector.size(); ++i) {
            m_QVector << qvector[i][0] <<","<< qvector[i][1] << " ";
        }
        m_QVector << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }

    if (m_HQVector.is_open()) {
        for (size_t i = 0; i < hqvector.size(); ++i) {
            m_HQVector << hqvector[i] << " ";
        }
        m_HQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
    if (m_HPQVector.is_open()) {
        for (size_t i = 0; i < hpqvector.size(); ++i) {
            m_HPQVector << hpqvector[i] << " ";
        }
        m_HPQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }

    if (m_PQVector.is_open()) {
        for (size_t i = 0; i < pqvector.size(); ++i) {
            m_PQVector << pqvector[i] << " ";
        }
        m_PQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }   
    }
    
}

void FluctuationSpectrum::CalculateSpectrumIndividualComplex(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    double Lx=(*(m_pState->GetMesh()->GetBox()))(0);
    double Ly=(*(m_pState->GetMesh()->GetBox()))(1);
    double Lz=(*(m_pState->GetMesh()->GetBox()))(2);

    m_Lz=Lz;
    AverageHeight();



    

    

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        const std::vector<vertex*>& pAllVertices = m_pState->GetMesh()->GetActiveV();
        int number_of_vertices = pAllVertices.size();
        int number_of_inclusions = (m_pState->GetMesh()->GetInclusion()).size(); // Assuming GetInclusions() returns the binclusion vector
        m_AverageInclusionDensity = static_cast<double>(number_of_inclusions) / number_of_vertices;
    }

    std::vector<std::vector<double>> qvector(m_MatrixOrderIndividual.size(), std::vector<double>(2, 0));
    std::vector<Complex> hqvector(m_MatrixOrderIndividual.size(),0);
    std::vector<Complex> pqvector(m_MatrixOrderIndividual.size(),0);

    for (size_t i = 0; i < m_MatrixOrderIndividual.size(); ++i) {
        std::vector<double> q={m_MatrixOrderIndividual[i][0]*2*PI/Lx,m_MatrixOrderIndividual[i][1]*2*PI/Ly};
        qvector[i]=q;
    }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
        FourierResult final_result=FourierTransformInclusionIndividualComplex(qvector);
        hqvector = final_result.sum_h;
        pqvector = final_result.sum_p;
    }
    else{
        std::vector<Complex> final_result=FourierTransformNoInclusionIndividualComplex(qvector);
        hqvector = final_result;
    }


    

    
    // Write the qvector and hqvector to files
    

    if (m_QVector.is_open()) {
        for (size_t i = 0; i < qvector.size(); ++i) {
            m_QVector << qvector[i][0] <<","<< qvector[i][1] << " ";
        }
        m_QVector << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }

    if (m_HQVector.is_open()) {
        for (size_t i = 0; i < hqvector.size(); ++i) {
            m_HQVector << hqvector[i].real() <<","<< hqvector[i].imag() << " ";
        }
        m_HQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }


    if (m_PQVector.is_open()) {
        for (size_t i = 0; i < pqvector.size(); ++i) {
            m_PQVector << pqvector[i].real() <<","<< pqvector[i].imag() << " ";
        }
        m_PQVector << "\n";
    } else {
        std::cerr << "Unable to open pqvector.txt for writing." << std::endl;
    }

}
    

void FluctuationSpectrum::CalculateSpectrumSphericalComplex(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    //1. Calculate Volume
    //2. Get center of mass
    double volume=m_pState->GetAnalysisCalculations()->GetVolume();
    double radius_0=std::cbrt((3.0*volume)/(4.0*PI));
    m_radius0=radius_0;
    CenterOfMass();

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()) {
        const std::vector<vertex*>& pAllVertices = m_pState->GetMesh()->GetActiveV();
        int number_of_vertices = pAllVertices.size();
        int number_of_inclusions = (m_pState->GetMesh()->GetInclusion()).size(); // Assuming GetInclusions() returns the binclusion vector
        m_AverageInclusionDensity = static_cast<double>(number_of_inclusions) / number_of_vertices;
    }

    std::vector<std::vector<int>> qvector(m_VectorOrder.size(), std::vector<int>(2, 0));
    std::vector<Complex> hqvector(m_VectorOrder.size(),0);
    std::vector<Complex> pqvector(m_VectorOrder.size(),0);

    for (size_t i = 0; i < m_VectorOrder.size(); ++i) {
        std::vector<int> q={m_VectorOrder[i][0],m_VectorOrder[i][1]};
        qvector[i]=q;
    }

    if (m_pState->GetAnalysisVariables()->GetInclusionCalculationActive()){
        FourierResult final_result=FourierTransformInclusionSphericalComplex(qvector);
        hqvector = final_result.sum_h;
        pqvector = final_result.sum_p;
    }
    else{
        FourierResult final_result=FourierTransformInclusionSphericalComplex(qvector);
        hqvector = final_result.sum_h;
        pqvector = final_result.sum_p;
    }


    

    
    // Write the qvector and hqvector to files
    

    if (m_QVector.is_open()) {
        for (size_t i = 0; i < qvector.size(); ++i) {
            m_QVector << qvector[i][0] <<","<< qvector[i][1] << " ";
        }
        m_QVector << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }

    if (m_HQVector.is_open()) {
        for (size_t i = 0; i < hqvector.size(); ++i) {
            m_HQVector << hqvector[i].real() <<","<< hqvector[i].imag() << " ";
        }
        m_HQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }


    if (m_PQVector.is_open()) {
        for (size_t i = 0; i < pqvector.size(); ++i) {
            m_PQVector << pqvector[i].real() <<","<< pqvector[i].imag() << " ";
        }
        m_PQVector << "\n";
    } else {
        std::cerr << "Unable to open pqvector.txt for writing." << std::endl;
    }

}


void FluctuationSpectrum::SphericalHarmonic(int l, int m, double theta, double phi, double& real_part, double& imag_part) {
    // Calculate associated Legendre polynomial P_l^m(cos(theta))
    double plm = AssociatedLegendre(l, std::abs(m), std::cos(theta));
    
    // Normalization factor
    double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * PI) * 
                           Factorial(l - std::abs(m)) / Factorial(l + std::abs(m)));
    
    // Calculate Y_l^m
    if (m == 0) {
        real_part = norm * plm;
        imag_part = 0.0;
    } else if (m > 0) {
        real_part = norm * plm * std::cos(m * phi) * std::sqrt(2.0);
        imag_part = norm * plm * std::sin(m * phi) * std::sqrt(2.0);
    } else { // m < 0
        real_part = norm * plm * std::sin(std::abs(m) * phi) * std::sqrt(2.0) * (std::abs(m) % 2 == 0 ? 1 : -1);
        imag_part = -norm * plm * std::cos(std::abs(m) * phi) * std::sqrt(2.0) * (std::abs(m) % 2 == 0 ? 1 : -1);
    }
}

// Add these helper functions before SphericalHarmonic

double FluctuationSpectrum::Factorial(int n) {
    if (n <= 1) return 1.0;
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

double FluctuationSpectrum::AssociatedLegendre(int l, int m, double x) {
    // Handle special cases
    if (m < 0 || m > l) return 0.0;
    if (l == 0) return 1.0;
    
    // Calculate P_l^m(x) using recurrence relations
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = std::sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    
    if (l == m) return pmm;
    
    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) return pmmp1;
    
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    
    return pll;
}