#include <iostream>
#include "Vec3D.h" 
#include <math.h>
#include <limits.h>
Vec3D::Vec3D(double x,double y,double z)
{
    m_X=x;
    m_Y=y;
    m_Z=z;
  
}
Vec3D::Vec3D()
{
    m_X=0.0;
    m_Y=0.0;
    m_Z=0.0;
}
Vec3D::~Vec3D()
{
   
}

double Vec3D::at(int n)
{
double Value=0;
    if(n==0)
    {
        Value=m_X;
    }
    else if(n==1)
    {
        Value=m_Y;
    }
    else if(n==2)
    {
        Value=m_Z;
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
    }
return Value; 
}
void Vec3D::put( int n, double s)
{

    if(n==0)
    {
        m_X=s;
    }
    else if(n==1)
    {
        m_Y=s;
    }
    else if(n==2)
    {
        m_Z=s;
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
    }
}
double& Vec3D::operator()(const int n)
{
    double *Value=0;
    if(n==0)
    {
        Value=&m_X;
    }
    else if(n==1)
    {
        Value=&m_Y;
    }
    else if(n==2)
    {
        Value=&m_Z;
    }
    else
    {
        std::cout<<"Error: index should not be larger the 2 \n";
    }
    return *Value;
}
//-------------------------------------------------------
Vec3D Vec3D::operator+(Vec3D M)
{
    Vec3D M1;
    M1(0)=this->at(0)+M(0);
    M1(1)=this->at(1)+M(1);
    M1(2)=this->at(2)+M(2);
    return M1;
}
//-------------------------------------------------------
Vec3D Vec3D::operator-(Vec3D M)
{
    Vec3D M1;
    M1(0)=this->at(0)-M(0);
    M1(1)=this->at(1)-M(1);
    M1(2)=this->at(2)-M(2);
    return M1;
}
//-------------------------------------------------------
Vec3D Vec3D::operator*(Vec3D M)
{
    Vec3D M1;
     M1(0)=(this->at(1))*M(2)-(this->at(2))*M(1);
     M1(1)=(this->at(2))*M(0)-(this->at(0))*M(2);
     M1(2)=(this->at(0))*M(1)-(this->at(1))*M(0);
    
    return M1;
}
Vec3D Vec3D::operator*(double x)
{
    Vec3D M1;
    M1(0)=(this->at(0))*x;
    M1(1)=(this->at(1))*x;
    M1(2)=(this->at(2))*x;

    return M1;
}
void Vec3D::operator=(Vec3D M)
{
            m_X= M(0);
            m_Y= M(1);
            m_Z= M(2);
}
double Vec3D::norm ()
{
    double No=m_X*m_X+m_Y*m_Y+m_Z*m_Z;
    return sqrt(No);
}
double Vec3D::dot(Vec3D v1,Vec3D v2)
{
    double No=v1(0)*v2(0)+v1(1)*v2(1)+v1(2)*v2(2);
    return No;
}
