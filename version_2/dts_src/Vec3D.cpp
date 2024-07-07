#include <iostream>
#include <math.h>
#include <limits.h>
#include "Vec3D.h" 
/*
 Perhaps we could use array so we do not have if condition for accessing membranes 
 */
std::ostream& operator<<(std::ostream& os, const Vec3D& vec) {
    os << vec.m_X << " " << vec.m_Y << " " << vec.m_Z;
    return os;
}
Vec3D::Vec3D(double x, double y, double z) : m_X(x), m_Y(y), m_Z(z) {}

double& Vec3D::operator()(const int n) {
   // if (n < 0 || n >= 3) throw std::out_of_range("Index out of range");
    
    if (n == 0) return m_X;
    else if (n == 1) return m_Y;
    else return m_Z;
}

Vec3D Vec3D::operator+(const Vec3D& other) const {
    return Vec3D(m_X + other.m_X, m_Y + other.m_Y, m_Z + other.m_Z);
}


Vec3D Vec3D::operator-(const Vec3D& other) const {
    return Vec3D(m_X - other.m_X, m_Y - other.m_Y, m_Z - other.m_Z);
}

Vec3D Vec3D::operator*(const Vec3D& other) const {
    return Vec3D(m_Y * other.m_Z - m_Z * other.m_Y,
                 m_Z * other.m_X - m_X * other.m_Z,
                 m_X * other.m_Y - m_Y * other.m_X);
}
void Vec3D::normalize(){
    double norm = std::sqrt(m_X*m_X + m_Y*m_Y + m_Z*m_Z);
    
    if(norm == 0){
        return;
    }
    
    m_X = m_X/norm;
    m_Y = m_Y/norm;
    m_Z = m_Z/norm;

    
}
Vec3D Vec3D::operator*(double scalar) const {
    return Vec3D(m_X * scalar, m_Y * scalar, m_Z * scalar);
}

void Vec3D::operator=(const Vec3D& other) {
    m_X = other.m_X;
    m_Y = other.m_Y;
    m_Z = other.m_Z;
}

double Vec3D::norm() const {
    return std::sqrt(m_X * m_X + m_Y * m_Y + m_Z * m_Z);
}
double Vec3D::dot(const Vec3D& v1, const Vec3D& v2) {
    return v1.m_X * v2.m_X + v1.m_Y * v2.m_Y + v1.m_Z * v2.m_Z;
}

bool Vec3D::isbad() const {
    
    return !std::isfinite(m_X) || !std::isfinite(m_Y) || !std::isfinite(m_Z);
}
bool Vec3D::isgood() const {
    return std::isfinite(m_X) && std::isfinite(m_Y) && std::isfinite(m_Z);
}
void Vec3D::print(){
    std::cout<<m_X<<"  "<<m_Y<<"  "<<m_Z<<"\n";
}

/*
 
 bool Vec3D::isbad() const {
     return std::isinf(m_X) || std::isnan(m_X) ||
            std::isinf(m_Y) || std::isnan(m_Y) ||
            std::isinf(m_Z) || std::isnan(m_Z);
 }
 
 */
