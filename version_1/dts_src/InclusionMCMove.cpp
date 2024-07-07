

#include <stdio.h>
#include "InclusionMCMove.h"
#include "Curvature.h"
#include "inclusion.h"
#include "State.h"
InclusionMCMove::InclusionMCMove ()
{
}
InclusionMCMove::InclusionMCMove (State *pState)  {
    m_pState = pState;
    m_pInt     = m_pState->m_pinc_ForceField;
    m_ptotenergy=&(m_pState->m_TotEnergy);
}
InclusionMCMove::~InclusionMCMove()
{
    
}
void InclusionMCMove::MC_Move_AnInclusion(inclusion *pinc, RNG* Random)
{
    m_Beta = m_pState->m_Beta;

    m_LIntEChange.clear();
    m_pLIntEChange.clear();
    
    double R=0.5;
    m_pInc = pinc;
    m_MoveValidity=0;

    double temp = Random->UniformRNG(1);
    double movetype = Random->UniformRNG(1);
    std::vector <links*> linklist=(m_pInc->Getvertex())->GetVLinkList();
    int m = Random->IntRNG(linklist.size());
    double dx = (1-2*(Random->UniformRNG(1)));
    double dy = (1-2*(Random->UniformRNG(1)));
    links *d_link = linklist.at(m);

    if(movetype>0.5)
        KawasakiMove(temp, d_link);
    else
       RotationMove(temp,R*dx,R*dy);
}
void InclusionMCMove::KawasakiMove(double temp, links * d_links)
{
    double m_oldEnergy=0.0;
    vertex * hver = m_pInc->Getvertex();   // vertex that own the inclsuion
    vertex *tver = d_links->GetV2();        // vertex that the inclusion may travel to
    bool hasit=tver->VertexOwnInclusion();   // check if the target vertex own any inclusion
    
#if TEST_MODE == Enabled
    if(hver->GetVID()==tver->GetVID())
        std::cout<<"  Error: KawasakiMove this should not happen \n";
#endif
    if(hasit==true && m_pInc->GetInclusionTypeID()==(tver->GetInclusion())->GetInclusionTypeID())
    {
        
        // if on the both vertex the inclsuion are same so do not perform any move
    }
    else
    {
        double v1benden=hver->GetEnergy();   // energy of the vertex with inclusion
        double v2benden=tver->GetEnergy();   // energy of the target vertex
        m_oldEnergy+=v1benden;              // add the energies to the old energy term
        m_oldEnergy+=v2benden;
        std::vector<links *> l1list=hver->GetVLinkList();     // links that their interaction energies may change
        std::vector<links *> l2list=tver->GetVLinkList();     // links that their interaction energies may change

        // We should clear this container
        m_pLIntEChange.insert(m_pLIntEChange.end(), l1list.begin(), l1list.end());   // a copy of the links
        
        // a copy of the links and avoiding to repeat the link between v1 and v2
                for (std::vector<links *>::iterator it = l2list.begin() ; it != l2list.end(); ++it)
                {
                    if((*it)->GetID()!=d_links->GetID() && ((*it)->GetMirrorLink())->GetID()!=d_links->GetID())
                        m_pLIntEChange.push_back(*it);
                }

#if TEST_MODE == Enabled
        int sizeL=m_pLIntEChange.size();
        for (int i=0;i<sizeL;i++)
        {
            for (int j=i+1;j<sizeL;j++)
            {
                if((m_pLIntEChange.at(i))->GetID()==(m_pLIntEChange.at(j))->GetID() || ((m_pLIntEChange.at(i))->GetMirrorLink())->GetID()==(m_pLIntEChange.at(j))->GetID())
                {
                    std::cout<<"Error1818: This should not happen, in the inclusion mc move, the container contains same objects several times \n";
                    std::cout<<(m_pLIntEChange.at(i))->GetID()<<"  "<<(m_pLIntEChange.at(j))->GetID()<<"  ";
                    std::cout<< ((m_pLIntEChange.at(i))->GetMirrorLink())->GetID()<<"  "<<(m_pLIntEChange.at(j))->GetID()<<"\n";
                }
                
            }
        }
#endif
        // add the interaction energies and make a copy to the links associated with them
        for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
        {
            m_oldEnergy+=2*((*it)->GetIntEnergy());
            m_LIntEChange.push_back(*(*it));
        }
        //perform the move
        if(hasit==true )
        {
            // give the host vertex the traget vertex inclusion (only is allowed if target vertex has one)
            (tver->GetInclusion())->Updatevertex(hver);
            hver->UpdateInclusion((tver->GetInclusion()));
            
        }
        else
        {
            hver->UpdateOwnInclusion(false);
        }
        // move the inclusion to the target vertex
        m_pInc->Updatevertex(tver);
        tver->UpdateInclusion(m_pInc);
        tver->UpdateOwnInclusion(true);
        
        
        Energy EE (m_pInt);
        
        double NewEnergy=0;
        
        
        // calculate the interaction energies
        for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
        {
            NewEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
        }
        
        // elastic energies
            NewEnergy+=EE.SingleVertexEnergy(hver);
            NewEnergy+=EE.SingleVertexEnergy(tver);
            m_EnergyDifference=NewEnergy-m_oldEnergy;
        // MC acceptance
        if(m_EnergyDifference<=0 )
        {
            
            (*m_ptotenergy)=(*m_ptotenergy)+m_EnergyDifference;
            m_MoveValidity=1;

        }
        else if(exp(-m_Beta*m_EnergyDifference)>temp)
        {
            (*m_ptotenergy)=(*m_ptotenergy)+m_EnergyDifference;
            m_MoveValidity=1;
        }
        else
        {
            // rejection
            
            if(hasit==true )
            {
              tver->UpdateInclusion(hver->GetInclusion());
              (tver->GetInclusion())->Updatevertex(tver);
            }
            else
                tver->UpdateOwnInclusion(false);
            
            
            m_pInc->Updatevertex(hver);
            hver->UpdateOwnInclusion(true);
            hver->UpdateInclusion(m_pInc);
            hver->UpdateEnergy(v1benden);
            tver->UpdateEnergy(v2benden);
            // give the original interaction energy to all the links and their mirrors 
            std::vector<links>::iterator itr=m_LIntEChange.begin();
            for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
            {
                double en=(*itr).GetIntEnergy();
                (*it)->UpdateIntEnergy(en);
                ((*it)->GetMirrorLink())->UpdateIntEnergy(en);
                
#if TEST_MODE == Enabled
                if((*itr).GetID()!=(*it)->GetID())
                    std::cout<<"  "<<(*itr).GetID()<<"  "<<(*it)->GetID()<<"error:KawasakiMove-  \n";
#endif
                ++itr;
                
            }
        }
        
    }
    
}
void InclusionMCMove::RotationMove(double temp, double dx, double dy)
{
    
{
        double m_oldEnergy = 0;
        vertex * ver = m_pInc->Getvertex();
        std::vector<links *> l_list=ver->GetVLinkList();
        m_pLIntEChange.insert(m_pLIntEChange.end(), l_list.begin(), l_list.end());

        for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
        {
            m_oldEnergy+=2*((*it)->GetIntEnergy());
            m_LIntEChange.push_back(*(*it));
        }
    
        double vbenden=ver->GetEnergy();
        m_oldEnergy+=vbenden;
    
        Vec3D LD = m_pInc->GetLDirection();
        Vec3D newLD(LD(0)+dx,LD(1)+dy,0);
        newLD=newLD*(1/newLD.norm());
        m_pInc->UpdateLocalDirection(newLD);
    
    Energy EE (m_pInt);
    double NewEnergy=0;
    
    
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
    {
        NewEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
    }
    
    
        NewEnergy+=EE.SingleVertexEnergy(ver);
    
    m_EnergyDifference=NewEnergy-m_oldEnergy;

    
    
    if(m_EnergyDifference<=0 )
    {
        
        (*m_ptotenergy)=(*m_ptotenergy)+m_EnergyDifference;
        m_MoveValidity=1;
        
    }
    else if(exp(-m_Beta*m_EnergyDifference)>temp)
    {
        (*m_ptotenergy)=(*m_ptotenergy)+m_EnergyDifference;
        m_MoveValidity=1;
    }
    else
    {
            m_pInc->UpdateLocalDirection(LD);
            ver->UpdateEnergy(vbenden);
            std::vector<links>::iterator itr=m_LIntEChange.begin();
            for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
            {
                    double en=(*itr).GetIntEnergy();
                    (*it)->UpdateIntEnergy(en);
                    ((*it)->GetMirrorLink())->UpdateIntEnergy(en);
            
#if TEST_MODE == Enabled
            
                    if((*itr).GetID()!=(*it)->GetID())
                        std::cout<<"  "<<(*itr).GetID()<<"  "<<(*it)->GetID()<<"error:Rotation move-  \n";
#endif
            ++itr;
            }
        
        
        
    }
    
}
}










