#include <iostream>
//#include <vector>
#include "LMatrix.h" 
#include <math.h>
#include <limits.h>

/*

3 5 2
2 1 4    M=2 N=3



*/
LMatrix::LMatrix(int M, int N)                // M row and N column 
{

m_M=M;
m_N=N;


//              initialize your LMatrix by 0
         for (int i=0;i<M;i++)
         {
          std::vector <double> R;
         for (int j=0;j<N;j++)
         {
         R.push_back(0.0);
         }
           
         m_Mat.push_back(R);
         }



  
}
LMatrix::LMatrix(int M, int N, char I)                // M row and N column
{
    m_Mat.clear();
    m_M=M;
    m_N=N;
    
    
    //              initialize your LMatrix by 0
    
    if (I=='I')
    for (int i=0;i<M;i++)
    {
        std::vector <double> R;
        for (int j=0;j<N;j++)
        {
            if(i==j)
            R.push_back(1.0);
            else
            R.push_back(0.0);
                
        }
        
        m_Mat.push_back(R);
    }
    
    
    
    
}



LMatrix::~LMatrix()
{
   
}







/////////////// Function

///////// F1



///////////////////////////  F3   x and y are the power of the x and y and last term is changing from 0 to  m_LMatrixTermno-1



double LMatrix::at( int n, int m)
{

double Value=0;


          std::vector <double> R=m_Mat.at(n);
          Value=R.at(m);


return Value; 
}
void LMatrix::put( int n, int m, double s)
{

  std::vector <double> R=m_Mat.at(n);
R=m_Mat.at(n);
R.at(m)=s;
m_Mat.at(n)=R;

}

void LMatrix::print( )
{




         for (int j=0;j<m_M;j++)
         {
         for (int i=0;i<m_N;i++)
         {
         std::cout<<((m_Mat.at(j)).at(i))<<" ";
         }
         
         std::cout<<"\n";
         }
 
}

double& LMatrix::operator()(const int n, const int m)
{


    return (m_Mat.at(n)).at(m);
          

}
//-------------------------------------------------------
LMatrix LMatrix::operator+(LMatrix M)
{
    LMatrix M1(m_M,m_N);
   // std::vector <std::vector <double> > m_Mat;


    if (this->Rsize()==M.Rsize() && this->Csize()==M.Csize())
    {
    for (int i=0;i<m_M;i++)
    {

        
        
        for (int j=0;j<m_N;j++)
        {
            M1(i,j)=this->at(i,j)+M(i,j);
        }
        
       
    }

   }
   else
   {
       std::cout<<" Error:  LMatrix are not same size, we can not add them then Sorry. No Action happend here "<<"\n";
    
   }
    return M1;
    
}
//-------------------------------------------------------
LMatrix LMatrix::operator-(LMatrix M)
{
    LMatrix M1(m_M,m_N);
    if (this->Rsize()==M.Rsize() && this->Csize()==M.Csize())
    {
    for (int i=0;i<m_M;i++)
    {

        
        
        for (int j=0;j<m_N;j++)
        {
            M1(i,j)=this->at(i,j)-M(i,j);
        }
        
       
    }

   }
   else
   {
       std::cout<<" Error:  LMatrix are not same size, we can not do substration.  Sorry. No Action happend here "<<"\n";
    
   }
    return M1;

}
//-------------------------------------------------------
LMatrix LMatrix::operator*(LMatrix M)
{

  if (this->Csize()==M.Rsize() )
    {
    int NN=M.Rsize();
    int rs=this->Rsize();
    int cs=M.Csize();
    LMatrix M1(rs,cs);
 

    for (int i=0;i<rs;i++)
    {
        
        
        
        for (int j=0;j<cs;j++)
        {
                for (int k=0;k<M.Rsize();k++)
                {

                     M1(i,j)=M1(i,j)+((this->at(i,k))*(M(k,j)));
                }
        }
        
       
    }
    
    return M1;
   }
   else
   {
       std::cout<<" Error:  LMatrix are not have condition to be multiplied.  Sorry. No Action happend here "<<"\n";
    return M;
   }
    
}
//-------------------------------------------------------
void LMatrix::operator=(LMatrix M)
{
   if (this->Rsize()==M.Rsize() && this->Csize()==M.Csize())
    {



    for (int i=0;i<m_M;i++)
    {
       
        
        
        for (int j=0;j<m_N;j++)
        {
            (m_Mat.at(i)).at(j)= M(i,j);
        }
        
        
    }

    }
      else
   {
       std::cout<<" Error:  LMatrix are not same size, they can not be equall "<<"\n";
    
   } 

    
}
//==========================================================
//=================== Inverse LMatrix =========================
//=============================================
LMatrix LMatrix::Inverse (LMatrix M)
{
    int rs=M.Rsize();
    int cs=M.Csize();
    LMatrix InM(rs,cs);
    if(rs!=cs)
    {
   std::cout<<"Error: there is no invearse for a LMatrix with Row=!coulmn. It returns an empty LMatrix "<<"\n";
    }
    else
    {
    LMatrix BM(3*rs,3*cs);

    
        for ( int i=0;i<rs;i++)
        {
            for ( int j=0;j<cs;j++)
            {
                BM(i+1,j+1)=M(i,j);
            }
            
        }
        
        

        
        for(int i=1;i<=rs;i++)
            for(int j=1;j<=2*rs;j++)
                if(j==(i+rs))
                    BM(i,j)=1;
        
        /************** partial pivoting **************/
        for(int i=rs;i>1;i--)
        {
            if(BM(i-1,1)<BM(i,1))
                for(int j=1;j<=rs*2;j++)
                {
                    double d=BM(i,j);
                    BM(i,j)=BM(i-1,j);
                    BM(i-1,j)=d;
                }
        }
        
        /********** reducing to diagonal  LMatrix ***********/
      
        for(int i=1;i<=rs;i++)
        {
            for(int j=1;j<=rs*2;j++)
                if(j!=i)
                {
                   double d=BM(j,i)/BM(i,i);
                    for(int k=1;k<=rs*2;k++)
                        BM(j,k)-=BM(i,k)*d;
                }
        }
        /************** reducing to unit LMatrix *************/
        for(int i=1;i<=rs;i++)
        {
            double d=BM(i,i);
            for(int j=1;j<=rs*2;j++)
                BM(i,j)=BM(i,j)/d;
        }
        
        
    
        
        for (int i=0;i<rs;i++)
        {
            
            
           
            for (int j=0;j<rs;j++)
            {
                InM(i,j)=BM(i+1,rs+1+j);
            }
            
            
        }
        

        
        
    }
        return InM;
 
    


}

LMatrix LMatrix::TM (LMatrix M)
{

int m = M.Rsize();
int n = M.Csize();
LMatrix T(n,m);


for (int i=0;i<m;i++)
{
for (int j=0;j<n;j++)
{

T(j,i)=M(i,j);

}


}





return T;
}
double LMatrix::norm ()
{
    
    
    double No=0;
    
    if(m_M == 1 )
    {
            for (int i=0;i<m_M;i++)
            {
                std::vector <double> R=m_Mat.at(i);
                for (int j=0;j<m_N;j++)
                {
                    No+=R.at(j);
                }
                
            }

    }
    else if(m_N == 1 )
    {
        for (int i=0;i<m_M;i++)
        {
            std::vector <double> R=m_Mat.at(i);
            for (int j=0;j<m_N;j++)
            {
                No+=R.at(j)*R.at(j);
            }
            
        }
        
    }
    else
    {
        std::cout<<"Error: the LMatrix is not a vector \n";
    }


    
    
    
    return sqrt(No);
}




//===========================================================================================
