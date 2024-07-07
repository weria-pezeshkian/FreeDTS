

#include <stdio.h>
#include "GenerateCNTCells.h"
#include "Nfunction.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Triangle object: points to 3 vertices and has an area and normal vector
 */
GenerateCNTCells::GenerateCNTCells(std::vector< vertex* > ver, Vec3D cntsize, Vec3D *pBox)
{
	m_pBox=pBox;
    m_Allver=ver;
    m_CNTSize.clear();
    m_CNTSize.push_back(cntsize(0));
    m_CNTSize.push_back(cntsize(1));
    m_CNTSize.push_back(cntsize(2));
}


GenerateCNTCells::~GenerateCNTCells()
{
    
}

int GenerateCNTCells::IDFromIndex(int i,int j,int k)
{

int n=i+m_Nx*m_Ny*k+m_Nx*j;

return n;
}

void GenerateCNTCells::Generate()
{
    
    m_pAllCNTCells.clear();
    m_AllCNTCells.clear();
    m_CNTCellSize.clear();
    m_CNTCellNo.clear();
    
    
    
    
    //==================
    //====
    //==================
    Nfunction f;
    std::vector <double> Box;
    Box.push_back((*m_pBox)(0));
    Box.push_back((*m_pBox)(1));
    Box.push_back((*m_pBox)(2));
    
    std::vector <double> CNTSize=m_CNTSize;
    
    
    
    /// The CNT cell should not be smaller then one; for dynamic box we should be carefull
    
    m_Nx=int(Box.at(0)/CNTSize.at(0));
    m_Ny=int(Box.at(1)/CNTSize.at(1));
    m_Nz=int(Box.at(2)/CNTSize.at(2));
    


    if(Box.at(0)/double(m_Nx)==CNTSize.at(0))
    {
        
    }
    else
    {
        //m_Nx=m_Nx+1; this is wrong it makes the box smaller 
        std::string A="CNT cell in x direction has the size of "+f.Int_to_String(Box.at(0)/double(m_Nx));
        std::cout<<A<<"\n";
        f.Write_One_LogMessage(A);
    }
    
    if(Box.at(1)/double(m_Ny)==CNTSize.at(1))
    {
        
    }
    else
    {
       // m_Ny=m_Ny+1;
        std::string A="CNT cell in y direction has the size of "+f.Int_to_String(Box.at(1)/double(m_Ny));
        std::cout<<A<<"\n";
        f.Write_One_LogMessage(A);
    }
    if(Box.at(2)/double(m_Nz)==CNTSize.at(2))
    {
        
    }
    else
    {
      //  m_Nz=m_Nz+1;
        std::string A="CNT cell in z direction has the size of "+f.Int_to_String(Box.at(2)/double(m_Nz));
        std::cout<<A<<"\n";
        f.Write_One_LogMessage(A);
    }
    
    
    
    
    m_CNTCellSize.push_back(Box.at(0)/double(m_Nx));
    m_CNTCellSize.push_back(Box.at(1)/double(m_Ny));
    m_CNTCellSize.push_back(Box.at(2)/double(m_Nz));
    m_CNTCellNo.push_back(m_Nx);
    m_CNTCellNo.push_back(m_Ny);
    m_CNTCellNo.push_back(m_Nz);
    
    
    
    
    
    int m=0;    // id of the CNT cell
    
    //======================================
    //====== Generate CNT
    //=======================================
    for (int i=0;i<m_Nz;i++)
    {
        for (int j=0;j<m_Ny;j++)
        {
            for (int k=0;k<m_Nx;k++)
            {
                CNTCell T(m,k,j,i);
                Vec3D Ss(k*(m_CNTCellSize.at(0)),j*(m_CNTCellSize.at(1)),i*(m_CNTCellSize.at(2)));
                Vec3D Sb((k+1)*(m_CNTCellSize.at(0)),(j+1)*(m_CNTCellSize.at(1)),(i+1)*(m_CNTCellSize.at(2)));
                T.UpdateCNTSidemin(Ss);
                T.UpdateCNTSidemax(Sb);
                m_AllCNTCells.push_back(T);
                m++;
            }
      	 }
    }
    
    
    //======================================
    //====== Adding Neighbour CNT Cell
    //=======================================
    m=0;
    for (int i=0;i<m_Nz;i++)
    {
        for (int j=0;j<m_Ny;j++)
        {
            for (int k=0;k<m_Nx;k++)
            {
                
                
                int id=0;
                
                int p=0;
                
                
                for (int sz=-1;sz<2;sz++)
                {
                    for (int sy=-1;sy<2;sy++)
                    {
                        for (int sx=-1;sx<2;sx++)
                        {
                            p++;
                            if(sz==0 && sx==0 &&sy==0)
                            {
                                
                            }
                            else
                            {
                                id=IDFromIndex((m_Nx+k+sx)%m_Nx,(m_Ny+j+sy)%m_Ny,(m_Nz+i+sz)%m_Nz);
                                (m_AllCNTCells.at(m)).AddtoNeighbourCNTCell(&(m_AllCNTCells.at(id)));
                            }
                        }
                    }
                }
                
                
                
                
                m++;
            }
            
        }
    }
    
    //======================================
    //====== Adding vertex to CNT
    //=======================================
    
    
    
    for (std::vector<vertex *>::iterator it = m_Allver.begin() ; it != m_Allver.end(); ++it)
    {
        
        double x=(*it)->GetVXPos();
        double y=(*it)->GetVYPos();
        double z=(*it)->GetVZPos();
        
        int nx=int(x/m_CNTCellSize.at(0));
        int ny=int(y/m_CNTCellSize.at(1));
        int nz=int(z/m_CNTCellSize.at(2));
        int id=IDFromIndex(nx,ny,nz);
        
        if(id>=m_AllCNTCells.size())
        {
            
            std::string A=" ERROR: particle with coordinate "+f.Int_to_String(x)+"  "+f.Int_to_String(y)+"  "+f.Int_to_String(z)+"  "+"   is out of the box \n";
            std::cout<<A<<"\n";
            f.Write_One_LogMessage(A);
            
        }
        else
        {
            (m_AllCNTCells.at(id)).AddtoVertexList(*it);
            (*it)->UpdateVCNTCell(&(m_AllCNTCells.at(id)));
        }
        
    }
    
    
    for (std::vector<CNTCell>::iterator it = m_AllCNTCells.begin() ; it != m_AllCNTCells.end(); ++it)
    {
        
        m_pAllCNTCells.push_back(&(*it));
    }
    
    
    
    
    
    
    
    
    

    
    
}



