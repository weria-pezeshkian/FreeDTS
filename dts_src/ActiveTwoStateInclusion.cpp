

#include <stdio.h>
#include "ActiveTwoStateInclusion.h"
#include "Nfunction.h"
#include "Energy.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Active inclusion
 This command is not tested and may contain errors, do not used it .....
 
 */

ActiveTwoStateInclusion::ActiveTwoStateInclusion()
{
    // ActiveTwoStateInclusion = on p1 p2 ep persentage gama
    m_Epsilon1 = 0;
    m_Epsilon2 = 0;
    m_Percentge = 0;
    m_State = false;
    m_Gama = 0;
    m_IncType1Name = "t1";
    m_IncType2Name = "t2";
    m_Health = true;
    m_N2 = 0;
    m_N1 = 0;
    m_N = 0;
    totDN = 0;
    m_ActiveEnergy = 0;;
    m_DeltaN = 0;
}
ActiveTwoStateInclusion::ActiveTwoStateInclusion(bool state, double ep1, double ep2,double per, double gama, std::string t1,std::string t2 )
{
    m_Epsilon1 = ep1;
    m_Epsilon2 = ep2;// at the moment this is not used anywhere...
    m_State = state;
    m_Percentge = per;
    m_Gama = gama;
    m_IncType1Name = t1;
    m_IncType2Name = t2;
    m_N2 = 0;
    m_N1 = 0;
    m_N = 0;
    totDN = 0;
    m_ActiveEnergy = 0;;
    m_DeltaN = 0;
    
    std::cout<<"Active membrane "<<m_IncType1Name<<"  "<<m_IncType2Name<<"  "<<m_Epsilon1<<" "<<m_Epsilon2<<"  "<<m_Percentge<<"  "<<m_Gama<<"\n";
}
void ActiveTwoStateInclusion::Initialize(std::vector<inclusion*> pInclusions, std::vector<InclusionType*> allinctype, Inclusion_Interaction_Map * pint, RNG* rng)
{
    
    /*
        m_Percentge is the difference between the two state; 1 and 2,
        m_N is total number of the inclusions that are involved in this exchanges. So DeltaN = 2*N1-m_N
        N1=int (m_Percentge*m_N)
        N2 = m_N - N1
        epsilon is the rate, we could have different rate but reducing the number of the model parameters 
     */
    Nfunction f;
    m_Health = true;
    m_RNG = rng;
    m_pInt = pint;
    m_pInclusions = pInclusions;
   if(m_State==true)
    {
        // adding all the relevant inclsuions into m_pInc
        for (std::vector<inclusion *>::iterator it = m_pInclusions.begin() ; it != m_pInclusions.end(); ++it)
        {
            if(((*it)->GetInclusionType())->ITName == m_IncType1Name)
            {
                m_N1++;
                m_pInc.push_back(*it);
            }
            else if(((*it)->GetInclusionType())->ITName == m_IncType2Name)
            {
                m_N2++;
                m_pInc.push_back(*it);
            }
        }
        // obtain a pointer to the inclusion type of the two inclusions
        for (std::vector<InclusionType *>::iterator it = allinctype.begin() ; it != allinctype.end(); ++it)
        {
            if((*it)->ITName == m_IncType1Name)
                m_pIncType1 = (*it);
            else if((*it)->ITName == m_IncType2Name)
                m_pIncType2 = (*it);
        }
    }
    m_N = m_pInc.size();            //total number of the involved inclusions
    int N10 = m_Percentge*m_N;      // averged number of the first type
    m_DeltaN0 = 2*N10-m_N;           // expected difference in the number of two states
    m_Gama  = m_Gama/m_N;           // fluctuations
    
    if (m_State==true)
    {
        std::string sms = " This system is coupled to an active process, i.e., ActiveTwoStateInclusion";
        f.Write_One_LogMessage(sms);
        std::cout<<sms<<std::endl;
        sms = " Active parameter are:\n inial N1 "+ f.Int_to_String(m_N1)+" inial N2 "+ f.Int_to_String(m_N2)+" target N_1-N2  "+ f.Int_to_String(m_DeltaN0);
        f.Write_One_LogMessage(sms);
        std::cout<<sms<<std::endl;

    }
}
ActiveTwoStateInclusion::~ActiveTwoStateInclusion()
{
    
}
void ActiveTwoStateInclusion::ActiveExchange(double *tot_Energy)
{
    
        Energy EE (m_pInt);
    
        double eta = 2*double(m_N-m_DeltaN0)/double(m_N+m_DeltaN0)-1; // we force ep1 to be equal to ep2. Does not make too much sense here
        // finding a random inclusion among the selected one

    for (int i=0;i<m_pInc.size();i++)
    {
        int  n = m_RNG->IntRNG(m_N);
        inclusion *tinc = m_pInc.at(n);
        vertex* tver = tinc->Getvertex();
        
        //==calculating the energies that may change
        std::vector<links *> Llist=tver->GetVLinkList();     // links that their interaction energies may change
        double oldEnergy=tver->GetEnergy(); // curvature energy of the vertex
        for (std::vector<links *>::iterator it = Llist.begin() ; it != Llist.end(); ++it)
            oldEnergy+=2*((*it)->GetIntEnergy()); // energy of the interacting inclusions (2 is for the mirror links which are not included)
        
        
        double newEnergy = oldEnergy;
        InclusionType* inctype = tinc->GetInclusionType();
        double thermal = m_RNG->UniformRNG(1.0);
        double phi1  = double(m_N1-m_N2-m_DeltaN0);
        double P = 0;
        if(inctype->ITName ==m_IncType1Name)// transition from 2-1
        {
            P =  m_Epsilon1*double(1)/double(m_N)*1.0/(1.0+exp(-m_Gama*phi1));// the difference between here and what we got in the paper is we loop around it and

            if(P>thermal)
            {
                tinc->UpdateInclusionType(m_pIncType2);
                tinc->UpdateInclusionTypeID(m_pIncType2->ITid);
                m_N1--;
                m_N2++;
                newEnergy = EE.SingleVertexEnergy(tinc->Getvertex());
                for (std::vector<links *>::iterator it = Llist.begin() ; it != Llist.end(); ++it)
                    newEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
            }

        }
        else  if(inctype->ITName == m_IncType2Name)// transition from 1-2
        {

            P =  m_Epsilon1*double(1)/double(m_N)*1.0/(eta+exp(m_Gama*phi1));
            if(P>thermal)
            {
                tinc->UpdateInclusionType(m_pIncType1);
                tinc->UpdateInclusionTypeID(m_pIncType1->ITid);
                m_N1++;
                m_N2--;
                newEnergy = EE.SingleVertexEnergy(tinc->Getvertex());
                for (std::vector<links *>::iterator it = Llist.begin() ; it != Llist.end(); ++it)
                    newEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
            }
        }
        (*tot_Energy)+= (newEnergy - oldEnergy);
        m_ActiveEnergy+= (newEnergy - oldEnergy);
    }
    

    m_DeltaN = m_N1-m_N2;

}
