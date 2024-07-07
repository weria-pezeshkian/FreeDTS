

#include <stdio.h>
#include "ActiveTwoStateInclusion.h"
#include "Nfunction.h"
#include "State.h"
 /*
  * @file ActiveTwoStateInclusion.cpp
  * @brief Implementation of the ActiveTwoStateInclusion class methods.
  *
  * This file contains the implementation of the methods declared in the ActiveTwoStateInclusion class.
  * The class is responsible for exchanging inclusion types between two states based on active algorthem.
  * The exchange is not done based on the energetic of the states. It is active, but energy get updated after the exchange
  * It initializes the inclusion exchange process, performs exchanges at defined intervals, and manages energy calculations.
  *
  * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
  */

ActiveTwoStateInclusion::ActiveTwoStateInclusion(int period, double ep, double per, double gama, std::string t_name1,std::string t_name2 ) :
                        m_Period(period),
                        m_Epsilon(ep),
                        m_Percentage(per),
                        m_Gama(gama),
                        m_TypeName_1(t_name1),
                        m_TypeName_2(t_name2),
                        m_N2(0),
                        m_N1(0)
{
    
}
void ActiveTwoStateInclusion::Initialize(State *pstate) {
    
    m_pState = pstate;
    const std::vector<inclusion *>& pAllInclusion = m_pState->GetMesh()->GetInclusion();
    for (std::vector<inclusion *>::const_iterator it = pAllInclusion.begin() ; it != pAllInclusion.end(); ++it){
        if((*it)->m_IncType->ITName == m_TypeName_1){
            m_pSubInc.push_back(*it);
            m_pIncType1 = (*it)->m_IncType;
            m_N1++;
        }
        else if((*it)->m_IncType->ITName == m_TypeName_2){
            m_pSubInc.push_back(*it);
            m_pIncType2 = (*it)->m_IncType;
            m_N2++;
        }
    } // end for (std::vector<inclusion *>::const_iterator it = pAllInclusion.begin() ; it != pAllInclusion.end(); ++it)
    
    m_Gama  = m_Gama/double(m_pSubInc.size());  // fluctuations, why do we divide it by N
    m_N = m_N1 + m_N2;
    m_ActiveEnergy = 0;
    m_Delta_N0 = m_N * (2 * m_Percentage - 1);   // expected difference in the number of two states; or forced
}
ActiveTwoStateInclusion::~ActiveTwoStateInclusion() {
    
}
std::string ActiveTwoStateInclusion::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}
bool ActiveTwoStateInclusion::Exchange(int step){
 
    
    if( step%m_Period !=0 )
        return false;
    
    // in the old version it was as: double eta = 2*double(m_N-m_DeltaN0)/double(m_N+m_DeltaN0)-1;
    double eta = 2*double(m_N2)/double(m_N1)-1; // we force ep1 to be equal to ep2. Does not make too much sense here

    for (std::vector<inclusion *>::iterator it = m_pSubInc.begin() ; it != m_pSubInc.end(); ++it){
        
        double old_energy = 0;
        double new_energy = 0;
        bool exchange_happend = false;
        int delta_N = m_N1 - m_N2;
        double phi = double(delta_N-m_Delta_N0);

        vertex* tver = (*it)->Getvertex();
        InclusionType* inctype = (*it)->m_IncType;
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        
        if((*it)->m_IncType == m_pIncType1){ // transition from 1->2
            double Prob = m_Epsilon / static_cast<double>(m_N) / (1.0 + exp(-m_Gama * phi));

            if( Prob > thermal ) {
                (*it)->m_IncType = m_pIncType2;
                m_N1--;
                m_N2++;
                exchange_happend = true;
            }
        } // if((*it)->m_IncType->ITName == m_TypeName_1)
        else { // transition from 2->1
            double Prob = m_Epsilon / static_cast<double>(m_N) / (eta + exp(m_Gama * phi));

            if(Prob > thermal) {
                (*it)->m_IncType  = m_pIncType1;
                m_N1++;
                m_N2--;
                exchange_happend = true;
            }
        }
        if(exchange_happend){ // note, while the exchange has happend, energies are not updated. so old energy can be obtained.
            
            old_energy += tver->GetEnergy();
            std::vector<links *> n_edges = tver->GetVLinkList();     // links that their interaction energies may change
            
            for (std::vector<links *>::iterator it = n_edges.begin() ; it != n_edges.end(); ++it){
                old_energy += 2*((*it)->GetIntEnergy());
            }
            //--- now, lets update energy
            new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(tver); //
            for (std::vector<links *>::iterator it = n_edges.begin() ; it != n_edges.end(); ++it){
                    new_energy += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            }
            //--- update the total elastic energy of the system
            double T_en = m_pState->GetEnergyCalculator()->GetEnergy();
            T_en += new_energy - old_energy;
            m_pState->GetEnergyCalculator()->UpdateTotalEnergy(T_en);
            m_ActiveEnergy += (new_energy - old_energy);
            
            m_NumberOfAttemptedMoves++;
            m_AcceptedMoves++;
            
        } // if(exchange_happend)
        else { // otherwise
            m_NumberOfAttemptedMoves++;
        }

        
    } //for (m_pSubInc.begin())
    
    
    return true;
}
