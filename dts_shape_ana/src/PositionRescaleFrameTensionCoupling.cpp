

#include <stdio.h>
#include "PositionRescaleFrameTensionCoupling.h"
#include "Nfunction.h"
#include "vertex.h"
#include "Curvature.h"
#include "State.h"

/*
=================================================================================================================
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "MCMoveBoxChange" function.
 
 This Class need to be optimised soon.............
=================================================================================================================
*/
PositionRescaleFrameTensionCoupling::PositionRescaleFrameTensionCoupling()
{
    
    
}
PositionRescaleFrameTensionCoupling::PositionRescaleFrameTensionCoupling(double sigmap,State *pstate)
{

    m_DetaR = 0;
    m_DeltaA = 0;
    m_pInc=pstate->m_pinc_ForceField;
    Nfunction f;
    m_pBox=(pstate->m_pMesh)->m_pBox;
    m_SigmaP=sigmap;

    m_dr=0.0;
    m_Lyx=(*m_pBox)(1)/(*m_pBox)(0);
    m_pLmin2 = &(pstate->m_MinVerticesDistanceSquare);
    m_pLmax2 = &(pstate->m_MaxLinkLengthSquare);
    m_pminAngle = &(pstate->m_MinFaceAngle);
    m_step = 0;
    m_pCFGC = pstate->GetGlobalCurvature();
    m_pSPBTG = pstate->Get2GroupHarmonic();
}


PositionRescaleFrameTensionCoupling::~PositionRescaleFrameTensionCoupling()
{
    
}
bool PositionRescaleFrameTensionCoupling::MCMoveBoxChange(double dr, double * tot_Energy, double temp, int step, GenerateCNTCells *pGenCNT, std::vector<vertex *> pAllVertex, std::vector<links *> pAllLinks, std::vector<triangle* > pAllTri)
{
m_pAllVertex = pAllVertex;
m_pAllLinks = pAllLinks;
m_pAllTriangle = pAllTri;
m_pGenCNT = pGenCNT;
m_step = step;

double nv=pAllVertex.size();

///////////   1. check the size of the cnt cell
///////////   2. check the distance  
///////////   3. if true copy the vertices, links and triangles 
///////////   4. update x,y
///////////   5. calulate area normal curvature and energy 
///////////   6. check the energy and accept or reject the move
    double DetaR = 0;
    double DeltaA = 0;
    


    if(m_pCFGC->GetState()==true)
    {

        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
        {
            std::vector<double> C=(*it)->GetCurvature();
            double area = (*it)->GetArea();
            DeltaA+=area;
            DetaR+=(C.at(0)+C.at(1))*area;

        }
    }

    
                    	m_Move=true;
    			        m_UpdateCNT=true;   // should be set to true. true means do not update it.

			CheckCNTSize();      // CNT cell should not be smaller then 1 after the move

			double DE=0;
			double DEArea=0;    // energy contribuion from area change 
            		m_dr=dr;
            		m_drx=dr;
            		m_dry=m_Lyx*dr;

           		m_oldLx=(*m_pBox)(0);
           		m_oldLy=(*m_pBox)(1);
           		m_newLx=(*m_pBox)(0)+m_drx;
           		m_newLy=(*m_pBox)(1)+m_dry;

    			m_Lnox = m_newLx/m_oldLx;
    			m_Lnoy = m_newLy/m_oldLy;
			double AreaRatio=(m_Lnox*m_Lnoy);
    
    


			if(m_Move==true)
    			{
        			
				CheckMinDistance();
    			}
			if(m_Move==true)
    			{
        			
				CheckLinkLength();
    			}
			if(m_Move==true)
    			{

			
        			PerformMove();

        			
			    for (std::vector<triangle *>::iterator it = m_pAllTriangle.begin() ; it != m_pAllTriangle.end(); ++it)
    				{

        				(*it)->UpdateNormal_Area(m_pBox);

    				}
			    for (std::vector<links *>::iterator it = m_pAllLinks.begin() ; it != m_pAllLinks.end(); ++it)
    				{

       					 (*it)->UpdateNormal();
       					 (*it)->UpdateShapeOperator(m_pBox);
  
    				}
   			    for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
    				{
      
        				Curvature P(*it);
                        
    				}

				    Energy En(m_pInc);
   
                    double newEnergy = En.TotalEnergy(m_pAllVertex, pAllLinks);
				    double oldEnergy = (*tot_Energy);
                    
                    double harmonicde=0;
                    if(m_pSPBTG->GetState()==true)
                    {
                        m_pSPBTG->CalculateEnergy(m_step);
                        double e1 = m_pSPBTG->GetEnergy();
                        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
                        {
                            double x = (*it)->GetVXPos();
                            double y = (*it)->GetVYPos();
                            Vec3D dr(m_Lnox*x-x,m_Lnoy*y-y,0);
                            m_pSPBTG->MovingVertex((*it), dr);
                        }
                        m_pSPBTG->CalculateEnergy(m_step);
                        double e2 = m_pSPBTG->GetEnergy();
                        harmonicde = e2-e1;
                    }

					       DE= newEnergy-oldEnergy+harmonicde;
                    
                //    std::cout<<newpullenergy-oldpullenergy<<" "<<oldpullenergy<<"\n";
                    
					 	DEArea = m_SigmaP*(m_newLx*m_newLy-m_oldLx*m_oldLy);

                    double eG = 0;
                    if(m_pCFGC->GetState()==true)
                    {
                        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
                        {
                            std::vector<double> C=(*it)->GetCurvature();
                            double area = (*it)->GetArea();
                            DeltaA-=area;
                            DetaR-=(C.at(0)+C.at(1))*area;
                        }

                        eG = m_pCFGC->CalculateEnergyChange(-m_DeltaA,-m_DetaR);

                    }
      						if(pow((AreaRatio),nv)*exp(-DE-DEArea-eG)>temp )
     						{
                                if(m_pCFGC->GetState()==true)
                                m_pCFGC->UpdateEnergyChange(-m_DeltaA,-m_DetaR);

      						}
      						else
      						{
          						RejectMove();
          						m_Move=false;
          
      						}
    			}
			if(m_Move==true)
    			{
        			
				CheckFaceAngle();
				if(m_Move==true)
    				{
        			
           						AcceptMove();
           						(*tot_Energy)=(*tot_Energy)+DE;
    				}
      				else
      				{
          						RejectMove();
          						m_Move=false;
          
      				}
    			}


			return m_Move;




    
}

/*
 ======================================================
 Checking for validity of the CNT Cells
 =========================================================
 */
void PositionRescaleFrameTensionCoupling::CheckCNTSize()
{

    m_pAllCNT = m_pGenCNT->GetAllCNTCells();
    Nfunction f;
    CNTCell*  TemCNT=m_pAllCNT.at(0);


    Vec3D cntlength=TemCNT->GetCNTSidemax()-TemCNT->GetCNTSidemin();

    if(cntlength(1)<1.8 || cntlength(0)<1.8)
    {
        m_UpdateCNT=false;
        m_Move=false;
        std::string sms=" CNT cell is small. C_Y= "+f.Int_to_String(cntlength(1))+" and C_X= "+f.Int_to_String(cntlength(0))+".  it is re-generated at step "+f.Int_to_String(m_step);
        f.Write_One_LogMessage(sms);
        
    }
    else
    {
        m_UpdateCNT=true;
    }

}


/*======================================================================
 This class goes through CNT cells  checks the distance between vertices after the box change. 
 ======================================================================*/
void PositionRescaleFrameTensionCoupling::CheckMinDistance()
{


    Vec3D Box;
    Box(0)=m_newLx;
    Box(1)=m_newLy;
    Box(2)=(*m_pBox)(2);
    
    for (std::vector<CNTCell *>::iterator it = m_pAllCNT.begin() ; it != m_pAllCNT.end(); ++it)
    {		


	if(m_Move==false)
		break;

	    std::vector <vertex *> pver=(*it)->GetVertexList();
	for (int i=0;i<pver.size();i++)
	{
		for (int j=i+1;j<pver.size();j++)
		{
		    double l2=DistanceSquardBetweenTwoVertices(pver.at(i),pver.at(j),Box);
		    if(l2<(*m_pLmin2))
                    {
                        m_Move=false;
                    }
		}

	}


            std::vector <CNTCell *> Nib=(*it)->GetVNeighbourCNTCell();
        for (std::vector<CNTCell *>::iterator it2 = Nib.begin() ; it2 != Nib.end(); ++it2)
        {
            
            std::vector <vertex *> pverN=(*it2)->GetVertexList();
            
            for (std::vector<vertex *>::iterator itv1 = pver.begin() ; itv1 != pver.end(); ++itv1)
            {
            for (std::vector<vertex *>::iterator itv2 = pverN.begin() ; itv2 != pverN.end(); ++itv2)
            {
                    
                    double l2=DistanceSquardBetweenTwoVertices(*itv1,*itv2,Box);
                    
                    if(l2<(*m_pLmin2))
                    {
                        m_Move=false;
                    }
            }
            }
            
            
            
            
        }




    }
    

    

    
}

double PositionRescaleFrameTensionCoupling::DistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2,Vec3D Box)
{
    
    

    
    double x2=v2->GetVXPos();
    double y2=v2->GetVYPos();
    double z2=v2->GetVZPos();
    double x1=v1->GetVXPos();
    double y1=v1->GetVYPos();
    double z1=v1->GetVZPos();
	x2=m_Lnox*x2;
	x1=m_Lnox*x1;
	y2=m_Lnoy*y2;
	y1=m_Lnoy*y1;


    double dx=x2-x1;
    double dy=y2-y1;
    double dz=z2-z1;
    
    
    if(fabs(dx)>Box(0)/2.0)
    {
        
        
        if(dx<0)
            dx=Box(0)+dx;
        else if(dx>0)
            dx=dx-Box(0);
        
        
    }
    if(fabs(dy)>Box(1)/2.0)
    {
        
        
        if(dy<0)
            dy=Box(1)+dy;
        else if(dy>0)
            dy=dy-Box(1);
        
        
    }
    
    if(fabs(dz)>Box(2)/2.0)
    {
        if(dz<0)
            dz=Box(2)+dz;
        else if(dz>0)
            dz=dz-Box(2);
    }
    
    
    double l2=dx*dx+dy*dy+dz*dz;
    


    
    return l2;
}
/*======================================================================
 This function goes through all links  checks the distance between vertices after the box change. 
 ======================================================================*/
void PositionRescaleFrameTensionCoupling::CheckLinkLength()
{
    Vec3D Box;
    Box(0)=m_newLx;
    Box(1)=m_newLy;
    Box(2)=(*m_pBox)(2);

    for (std::vector<links *>::iterator it = m_pAllLinks.begin() ; it != m_pAllLinks.end(); ++it)
    {

		vertex *v1= (*it)->GetV1();
		vertex *v2= (*it)->GetV2();

                    double l2=DistanceSquardBetweenTwoVertices(v1,v2,Box);
                    
                    if(l2>(*m_pLmax2))
                    {
                        m_Move=false;
                    }
  
    }

}
/*======================================================================
function to check the angle between two faces of a link
 ======================================================================*/
void PositionRescaleFrameTensionCoupling::CheckFaceAngle()
{

			    for (std::vector<links *>::iterator it = m_pAllLinks.begin() ; it != m_pAllLinks.end(); ++it)
    				{

		       				if (CheckFaceAngle((*it))==false)
						{
							m_Move=false;
							break;
						}
    				}

}

/*======================================================================
 This function makes a copy of vertices, links and trinagles incase of rejection. 
 ======================================================================*/
void PositionRescaleFrameTensionCoupling::PerformMove()
{

m_AllVertex.clear();
m_AllLinks.clear();
m_AllProjectedLinks.clear();
m_AllTriangle.clear();

(*m_pBox)(0)=m_newLx;
(*m_pBox)(1)=m_newLy;


        for (std::vector<links *>::iterator it = m_pAllLinks.begin() ; it != m_pAllLinks.end(); ++it)
        {

		m_AllLinks.push_back(*(*it));

		links* pro= (*it)->GetMirrorLink();
		m_AllProjectedLinks.push_back(*pro);

        }


        for (std::vector<triangle *>::iterator it = m_pAllTriangle.begin() ; it != m_pAllTriangle.end(); ++it)
        {

		m_AllTriangle.push_back(*(*it));

        }

        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
        {

		m_AllVertex.push_back(*(*it));

       }


        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
        {
		  double x = (*it)->GetVXPos();  
		  (*it)->UpdateVXPos(m_Lnox*x); 
		  double y = (*it)->GetVYPos();  
		  (*it)->UpdateVYPos(m_Lnoy*y);


        }



}



/*
 ======================================================
Function to accept the move and update everything
 =========================================================
 */
void PositionRescaleFrameTensionCoupling::AcceptMove()
{
    for (std::vector<CNTCell *>::iterator it = m_pAllCNT.begin() ; it != m_pAllCNT.end(); ++it)
    {
        Vec3D side1=(*it)->GetCNTSidemax();
        Vec3D side2=(*it)->GetCNTSidemin();
        side1(0)=side1(0)*m_Lnox;
        side1(1)=side1(1)*m_Lnoy;

        side2(0)=side2(0)*m_Lnox;
        side2(1)=side2(1)*m_Lnoy;


         (*it)->UpdateCNTSidemax(side1);
         (*it)->UpdateCNTSidemin(side2);
        
    }

}

void PositionRescaleFrameTensionCoupling::RejectMove()
{

    if(m_pSPBTG->GetState()==true)
    {
        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
        {
            double X = (*it)->GetVXPos();
            double Y = (*it)->GetVYPos();
            Vec3D dr(X-X/m_Lnox,Y-Y/m_Lnoy,0);
            m_pSPBTG->RejectMovingVertex((*it), dr);
        }
        m_pSPBTG->CalculateEnergy(m_step);
    }
    
    
    (*m_pBox)(0)=(*m_pBox)(0)-m_drx;
    (*m_pBox)(1)=(*m_pBox)(1)-m_dry;
	

	int i=0;
		  double x = (m_pAllVertex.at(0))->GetVXPos();  
        for (std::vector<vertex *>::iterator it = m_pAllVertex.begin() ; it != m_pAllVertex.end(); ++it)
        {
                *(*it)=m_AllVertex.at(i);
                i++;
        }

		  double x1 = (m_pAllVertex.at(0))->GetVXPos();  

        i=0;
        for (std::vector<links *>::iterator it = m_pAllLinks.begin() ; it != m_pAllLinks.end(); ++it)
        {
            *(*it)=m_AllLinks.at(i);
            links *pro=(*it)->GetMirrorLink();
			*pro = m_AllProjectedLinks.at(i);
            i++;
        }


	i=0;

        for (std::vector<triangle *>::iterator it = m_pAllTriangle.begin() ; it != m_pAllTriangle.end(); ++it)
        {
		*(*it)=m_AllTriangle.at(i);  
		 i++;
        }
}

bool   PositionRescaleFrameTensionCoupling::CheckFaceAngle(links * l)
{
    bool is=true;
    
   
    double t=0;
    
    Vec3D n=(l->GetTriangle())->GetNormalVector();
    Vec3D n3=((l->GetMirrorLink())->GetTriangle())->GetNormalVector();
    
    
    if(n.dot(n,n3)<(*m_pminAngle))
    {
        is=false;
    }

    
    
    return is;
}













