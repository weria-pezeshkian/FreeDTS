#include <iostream>
//#include <vector>
#include "Tensor2.h" 
#include <math.h>
#include <limits.h>

/*
    | V1 |
 T= | V2 |
    | V3 |
*/


Tensor2::Tensor2(Vec3D x,Vec3D y,Vec3D z)
{
    m_V1=x;
    m_V2=y;
    m_V3=z;
}
Tensor2::Tensor2()
{
}
Tensor2::Tensor2(char t)                // M row and N column
{
    if(t=='I')
    {
        m_V1(0)=1.0;
        m_V2(1)=1.0;
        m_V3(2)=1.0;
    }
    else if(t=='O')
    {
        m_V1(0)=0.0;
        m_V2(0)=0.0;
        m_V3(0)=0.0;
        m_V1(1)=0.0;
        m_V2(1)=0.0;
        m_V3(1)=0.0;
        m_V1(2)=0.0;
        m_V2(2)=0.0;
        m_V3(2)=0.0;

    }
    else
    {
        std::cout<<" Error: Matrix with name "<<t<<" has not been defined \n";
    }
}
Tensor2::~Tensor2()
{
   
}
double Tensor2::at(int n,int m)  // n row index m coulmn index
{

double Value=0;
    if(n==0)
    {
        Value=m_V1(m);
    }
    else if(n==1)
    {
        Value=m_V2(m);
    }
    else if(n==2)
    {
        Value=m_V3(m);
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
    }
return Value; 
}
void Tensor2::put( int n, int m, double s)
{
    if(n==0)
    {
        m_V1(m)=s;
    }
    else if(n==1)
    {
        m_V2(m)=s;
    }
    else if(n==2)
    {
        m_V3(m)=s;
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
    }
}
Vec3D Tensor2::MULT(Tensor2 T,Vec3D A)
{
    Vec3D v1=T.GetV1();
    Vec3D v2=T.GetV2();
    Vec3D v3=T.GetV3();
    double x=A.dot(v1,A);
    double y=A.dot(v2,A);
    double z=A.dot(v3,A);
    Vec3D V(x,y,z);
    
    return V;
}
Vec3D Tensor2::operator*(Vec3D A)
{
    double x=A.dot(m_V1,A);
    double y=A.dot(m_V2,A);
    double z=A.dot(m_V3,A);
    Vec3D V(x,y,z);
    return V;
}
double& Tensor2::operator()(const int n, const int m)
{

    double *Value=0;
    if(n==0)
    {
        Value=&(m_V1(m));
    }
    else if(n==1)
    {
        Value=&m_V2(m);
    }
    else if(n==2)
    {
        Value=&m_V3(m);
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
       
    }

    return *Value;
    
}
//-------------------------------------------------------
Tensor2 Tensor2::operator+(Tensor2 M)
{
    Tensor2 M1;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            M1(i,j)=this->at(i,j)+M(i,j);

        }
    }
    return M1;
}
Tensor2 Tensor2::operator-(Tensor2 M)
{
    Tensor2 M1;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            M1(i,j)=this->at(i,j)-M(i,j);
            
        }
    }
    return M1;
}
//-------------------------------------------------------
Tensor2 Tensor2::operator*(Tensor2 M)
{
    Tensor2 M1;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
            M1(i,j)=M1(i,j)+(this->at(i,k))*M(k,j);
            }
            
        }
    }
    return M1;
}
Tensor2 Tensor2::operator*(double x)
{
    Vec3D   v1=m_V1*x;
    Vec3D   v2=m_V2*x;
    Vec3D   v3=m_V3*x;
    Tensor2 M(v1,v2,v3);
    return M;
}

//-------------------------------------------------------
void Tensor2::operator=(Tensor2 M)
{
    for (int i=0;i<3;i++)
    {
        m_V1(i)=M(0,i);
         m_V2(i)=M(1,i);
         m_V3(i)=M(2,i);
    }
}
Tensor2 Tensor2::makeTen(Vec3D X)
{
    Tensor2 A;
    for (int i=0;i<3;i++)
    {
        A(0,i)=X(0)*X(i);
        A(1,i)=X(1)*X(i);
        A(2,i)=X(2)*X(i);
    }
    return A;
}

Tensor2 Tensor2::Transpose(Tensor2 X)
{
        Tensor2 A;
        A(0,1)=X(1,0);
        A(0,2)=X(2,0);
        A(1,0)=X(0,1);
        A(1,2)=X(2,1);
        A(2,0)=X(0,2);
        A(2,1)=X(1,2);
        A(0,0)=X(0,0);
        A(1,1)=X(1,1);
        A(2,2)=X(2,2);
        return A;
}
