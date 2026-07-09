

#include <stdio.h>
#include "ParallelTemperingMoveSimple.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#ifdef MPI_DETECTED
#include <mpi.h>

#include <sstream>



/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
ParallelTemperingMoveSimple::ParallelTemperingMoveSimple(State *pState,int period, int nprocessors, double minbeta, double maxbeta,int initialsteps) :
        m_pState(pState),
        m_Period(period),
        m_Nprocessors (nprocessors),
        m_MinBeta (minbeta),
        m_MaxBeta (maxbeta),
        m_InitialSteps(initialsteps)
        {
        }
            //-- convert the type string into the direction vector

ParallelTemperingMoveSimple::~ParallelTemperingMoveSimple() {
    
}
void ParallelTemperingMoveSimple::Initialize() {
    
    std::cout<<"---> the algorithm for Parallel Tempering involves applying this: "<< GetBaseDefaultReadName()<<" \n";
    #ifdef MPI_DETECTED
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"Rank: "<<rank<<std::endl;
    m_Rank=rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout<<"Size: "<<size<<std::endl;
    m_Size=size;
    m_TempID=rank;
    m_Counter=0;
    m_SizeIsEven=(m_Size%2==0);
    m_ReceiveBroadcast.resize(m_Size);
    m_RequestBroadcast.resize(m_Size);
    //m_RankWithUpTempID=rank;
    //m_RankWithDownTempID=rank;
    for (int i=0;i<m_Size;i++)
        {
            m_RankAtTempID.push_back(i);
        }
    
    std::vector<double> localBetaVec;
    if (rank == 0) {
        localBetaVec = ReadTemperatures();
    }

    // Broadcast the size of the vector to all ranks
    int vecSize = (rank == 0) ? localBetaVec.size() : 0;
    MPI_Bcast(&vecSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the local vector based on the size received
    std::vector<double> globalBetaVec(vecSize);
    if (rank == 0) {
        // Ensure localBetaVec.data() is the correct pointer for broadcasting
        MPI_Bcast(localBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(globalBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Now each rank has a copy of the vector
    m_BetaVec = globalBetaVec;
    if (rank==0){m_BetaVec=localBetaVec;}

    if (m_Size !=m_BetaVec.size()){
        std::cout<<"The selected temperatures do not match the number of processors."<<std::endl;
        exit(0);
    }


    m_pState->GetSimulation()->SetBeta(m_BetaVec[m_TempID], 0);
    //--> set the run tag id, we need to update this ID each time that the processor changes its temprature. The id should be temprature dependent
        //std::string gfile = ReplicaState.GetRunTag() + Nfunction::Int_to_String(beta); // general output file name
        //ReplicaState.UpdateRunTag(gfile);

    //In here, I need to think how to implement the restart in the right way (will be hard but it is extremly important to do it right)
    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(m_TempID));
    m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_TempID)+TimeSeriDataExt);
    m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(m_TempID));
    
    if (m_Restart==false){
        if (m_Rank == 0) {
            m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::out);
            m_TimeSeriesFile << "Rank-Temperature ID Mapping: "<<std::endl;
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile
        }
    }
    else{
        if (m_Rank == 0) {
            m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::app);
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile
        }

    }

    #endif
}

void ParallelTemperingMoveSimple::WriteRankAtTempIDToFile() {
    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "Error: Time series file is not open!" << std::endl;
        return;
    }

    for (const auto& rankID : m_RankAtTempID) {
        m_TimeSeriesFile << rankID << " ";
    }
    m_TimeSeriesFile << std::endl;
}


bool ParallelTemperingMoveSimple::EvolveOneStep(int step){
    /**
     * @brief a call function to exchange temperature between replicas.
     *
     * This function attempts a temperature swap between different replicas.
     * The algorithm is a decentralized one, in the sense that each couple of
     * neighbouring replicas that attempt the exchange do not need of the other ones
     * to continue with the simulation.
     * 
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
//---> if does not match the preiod, return false
    if(step%m_Period != 0 || step<m_InitialSteps)
        return false;

    
    bool CountIsEven=(m_Counter%2==0);

    //TEST::
    //CountIsEven=true;


    //RECEIVE BROADCAST FROM LAST EXCHANGE, SWAP RANKS, UPDATE RANKATTEMPID
    if (m_Counter>0){

        if (!CountIsEven){

            for(int i=0; i<m_Size; i++){
                if(i%2==0){
                    MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                }
            }

            for(int i=0; i<m_Size; i++){
                if(i%2==0){
                    if (m_ReceiveBroadcast[i]== 1){
                        std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                    }
                }
            }
            
            for(int i=0; i<m_Size; i++){
                if(m_RankAtTempID[i]==m_Rank){
                    m_TempID=i;
                    if (m_TempID==0){m_TargetState=true;}
                    break;
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", RankAtTempID:";
            //for (const int& element : m_RankAtTempID) {
            //    std::cout << element << ' ';
            //}
            //std::cout << std::endl; 

        }

        else if (CountIsEven){

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<"Checkpoint 1"<<std::endl;     

            for(int i=0; i<m_Size; i++){
                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                    }
                }

                else if (!m_SizeIsEven){
                    if(i%1==0){
                        MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                    }
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<"Checkpoint 2"<<std::endl;
            

            for(int i=0; i<m_Size; i++){

                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        if (m_ReceiveBroadcast[i]== 1){
                            std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                        }
                    }
                }
                else if (!m_SizeIsEven){
                    if(i%1==0){
                        if (m_ReceiveBroadcast[i]== 1){
                            std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                        }
                    }

                }
            }
            
            for(int i=0; i<m_Size; i++){
                if(m_RankAtTempID[i]==m_Rank){
                    m_TempID=i;
                    if (m_TempID==0){m_TargetState=true;}
                    break;
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", RankAtTempID:";
            //for (const int& element : m_RankAtTempID) {
            //    std::cout << element << ' ';
            //}
            //std::cout << std::endl; 

        }
        
    }
    //Just something
    std::fill(m_ReceiveBroadcast.begin(), m_ReceiveBroadcast.end(), 0);

    std::fill(m_RequestBroadcast.begin(), m_RequestBroadcast.end(), MPI_REQUEST_NULL);

    bool TempIDIsEven=(m_TempID%2==0);

    if (m_Rank==0){
        WriteRankAtTempIDToFile();
    }

    //ATTEMPT EXCHANGE
    //CountIsEven=true;
    if (CountIsEven){

    //IF COUNT IS EVEN, ODD TEMPIDS SEND ENERGY AND RECEIVE RESULT OF ATTEMPTED EXCHANGE
    if (!TempIDIsEven){

        //CALCULATE AND SEND ENERGY
        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        int ReceiveExchangeAcceptedInt;
        MPI_Send(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD);

        //RECEIVE RESULT OF EXCHANGE ATTEMPTS
        MPI_Recv(&ReceiveExchangeAcceptedInt, 1, MPI_INT, DestRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bool ReceiveExchangeAccepted=(ReceiveExchangeAcceptedInt==1);

        //IF ACCEPTED, EFFECTIVELY SWAP RANKS
        if (ReceiveExchangeAccepted){
            int NewTempID=DestTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID)); 
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt=0;
            MPI_Send(&EmptyBlockingInt, 1, MPI_INT, DestRank, 3, MPI_COMM_WORLD);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
            }
            
        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        for(int i=0;i<m_Size;i++){
            if(i%2==0){
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }


    }
    //IF COUNT IS EVEN, EVEN TEMPIDS RECEIVE ENERGY AND ATTEMPT EXCHANGE
    else if (TempIDIsEven){
        //DEALING WITH SPECIAL CASE. IF SIZE IS ODD, THE LAST TEMP ID IS OUT OF THE ATTEMPT.
        if (!m_SizeIsEven && m_TempID==m_Size-1){
            for(int i=0;i<m_Size;i++){
                if(i%2==0){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }

            m_Counter++;

            return false;
        }


        //CALCULATING OWN ENERGY AND RECEIVING ENERGY FROM HIGHER CLOSER ODD RANK (TEMPID +1)
        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Recv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        //ATTEMPT EXCHANGE
        double dE=receive_energy-own_energy;
        double dBeta=m_BetaVec[SourceTempID] - m_BetaVec[m_TempID];
        double RandomNumber=m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        //std::cout<<RandomNumber<<std::endl;
        bool ExchangeAccepted=(RandomNumber<exp(dBeta*dE));


        ExchangeAccepted=true;
        int ExchangeAcceptedInt = ExchangeAccepted ? 1 : 0;

        //SENDING RESULT OF ATTEMPT
        MPI_Send(&ExchangeAcceptedInt, 1, MPI_INT, SourceRank, 1, MPI_COMM_WORLD);
        
        
        //EFFECTIVE SWAP OF TEMPERATURES
        if (ExchangeAccepted){
            int NewTempID=SourceTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID));
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt;
            MPI_Recv(&EmptyBlockingInt, 1, MPI_INT, SourceRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
        }
        

        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        for(int i=0;i<m_Size;i++){
            if(i%2==0 && i!=m_TempID){
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
            else if(i==m_TempID){
                MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }

        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

    }

    }
    //IF COUNT IS ODD, ODD TEMPIDS RECEIVE ENERGY FROM TEMPID +1 AND ATTEMPT EXCHANGE.
    else if(!CountIsEven){

    if (!TempIDIsEven){
        //SPECIAL CASE: IF SIZE IS EVEN AND COUNT IS ODD, THE LAST TEMP ID IS OUT OF THE ATTEMPT.
        if (m_SizeIsEven && m_TempID==m_Size-1){
            for(int i=0;i<m_Size;i++){
                if(i%1==0 && i!=m_TempID){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            m_Counter++;
            //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step."<<std::endl;
            return false;
        }

        //CALCULATING OWN ENERGY AND RECEIVING ENERGY FROM HIGHER CLOSER ODD RANK (TEMPID +1)
        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Recv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        //ATTEMPT EXCHANGE
        double dE=receive_energy-own_energy;
        double dBeta=m_BetaVec[SourceTempID] - m_BetaVec[m_TempID];
        double RandomNumber=m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        //std::cout<<RandomNumber<<std::endl;
        bool ExchangeAccepted=(RandomNumber<exp(dBeta*dE));


        ExchangeAccepted=true;
        int ExchangeAcceptedInt = ExchangeAccepted ? 1 : 0;

        //SENDING RESULT OF ATTEMPT
        MPI_Send(&ExchangeAcceptedInt, 1, MPI_INT, SourceRank, 1, MPI_COMM_WORLD);
        
        
        //EFFECTIVE SWAP OF TEMPERATURES
        if (ExchangeAccepted){
            int NewTempID=SourceTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID));
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt;
            MPI_Recv(&EmptyBlockingInt, 1, MPI_INT, SourceRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
        }
        

        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        for(int i=0;i<m_Size;i++){
            if (!m_SizeIsEven){
                if(i%1==0 && i!=m_TempID ){

                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
                else if(i==m_TempID){
                    MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            else if(m_SizeIsEven){
                if(i%1==0 && i!=m_TempID && i!=m_Size-1){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
                else if(i==m_TempID){
                    MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
        }

        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

        
    }
    //ODD COUNT, EVEN TEMPID SEND ENERGY TO TEMPID -1.
    else if (TempIDIsEven){
        if (m_TempID==0){
            for(int i=0;i<m_Size;i++){
                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                    }
                }
                else if (!m_SizeIsEven){
                    if(i%1==0){
                        MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                    }
                }
            }
            //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step"<<std::endl;
            m_Counter++;
            return false;
        }


        //CALCULATE AND SEND ENERGY
        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        int ReceiveExchangeAcceptedInt;
        
        MPI_Send(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD);

        //RECEIVE RESULT OF EXCHANGE ATTEMPTS
        MPI_Recv(&ReceiveExchangeAcceptedInt, 1, MPI_INT, DestRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bool ReceiveExchangeAccepted=(ReceiveExchangeAcceptedInt==1);

        //IF ACCEPTED, EFFECTIVELY SWAP RANKS
        if (ReceiveExchangeAccepted){
            int NewTempID=DestTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID)); 
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt=0;
            MPI_Send(&EmptyBlockingInt, 1, MPI_INT, DestRank, 3, MPI_COMM_WORLD);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
            }
            
        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        for(int i=0;i<m_Size;i++){
            if (!m_SizeIsEven){
                if(i%1==0){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            else if(m_SizeIsEven){
                if(i%1==0 && i!=m_Size-1){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
        }
        
    }


    }
    //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step"<<std::endl;
    m_Counter++;
    return true;
}

//double own_energy=m_pState->GetEnergyCalculator()->GetEnergy(); **1
//double other_energy;
//std::cout<<"Energy Correctly Obtained: "<<own_energy<<std::endl;
//m_pState->GetRandomNumberGenerator()->UniformRNG(1.0) **2

bool ParallelTemperingMoveSimple::ChangeToNewTemperatureID(int NewTempID){
    //All that one must do when it is time to change Temperature ID
    //Set New Temp ID
    
    m_TempID=NewTempID;
    //Set New Beta
    m_pState->GetSimulation()->SetBeta(m_BetaVec[m_TempID], 0);
    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetFolderName() +"_" + Nfunction::Int_to_String(m_TempID));

    //Here, I close the file I was writing to
    //I move to the file I want to write to

    //I have to look into this. Do I have to close the file and open it again?
    //I haveto close the file and open it again...

    //In here I have to look at... what happens if one rank closes one file and opens another.
    //
    //m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_TempID)+TimeSeriDataExt);
    



    return true;
}

std::vector<double> ParallelTemperingMoveSimple::ReadTemperatures(){

    Nfunction f;
    std::string filename=ParallelTemperingMoveSimple::GetInitialTemperaturesFileName();

    if(f.FileExist(filename)==false)
    {
        std::cout<<"-----> Error: the temperature file name with the name "<<filename<<" does not exist"<<std::endl;
        exit(0);
    }

    std::ifstream TemperatureFile(filename);
    
    // Check if the file was opened successfully
    if (!TemperatureFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(0);
    }

    std::vector<double> betavec;

    // File processing code goes here

    std::string line;
    while (std::getline(TemperatureFile, line)) {
        // Skip empty or whitespace-only lines
        if (is_line_empty(line)) {
            continue;
        }

        std::stringstream ss(line);  // Use a stringstream to parse the line
        double value;
        // Read the double value from the line
        if (ss >> value) {
            // Ensure there is no extra data in the line
            if (ss >> std::ws && !ss.eof()) {
                std::cerr << "Error: More than one value or extra data in line: \"" << line << "\"" << std::endl;
                exit(0);  // Exit with error code
            }

            
            betavec.push_back(value);  // Add each double to the data vector
        } else {
            std::cerr << "Error: Could not read a double value from line: \"" << line << "\"" << std::endl;
            exit(0);  // Exit with error code
        }
    }

    // Close the file
    TemperatureFile.close();

    for (size_t i = 1; i < betavec.size(); ++i) {
        if (betavec[i] < betavec[i - 1]) {
            std::cerr << "Error: Vector is not sorted in ascending order." << std::endl;
            exit(0);  // Exit with error code
        }
    }

    for (size_t i = 0; i < betavec.size(); ++i) {
        betavec[i]=1/betavec[i];
        std::cout<<betavec[i]<<std::endl;
    }

    return betavec;

}

bool ParallelTemperingMoveSimple::is_line_empty(const std::string& line) {
    // Check if the line is empty or contains only whitespace
    return line.empty() || std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); });
}

bool ParallelTemperingMoveSimple::GetTargetState(){
    return m_TargetState;
}

void ParallelTemperingMoveSimple::SetRestart(){
    m_Restart=true;
}

std::string ParallelTemperingMoveSimple::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}


#endif