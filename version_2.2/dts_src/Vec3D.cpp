#include "Vec3D.h"

// Constructor
Vec3D::Vec3D(double x, double y, double z) : m_data{x, y, z} {}

// Addition
Vec3D Vec3D::operator+(const Vec3D& other) const {
    return Vec3D(m_data[0] + other.m_data[0], m_data[1] + other.m_data[1], m_data[2] + other.m_data[2]);
}

// Subtraction
Vec3D Vec3D::operator-(const Vec3D& other) const {
    return Vec3D(m_data[0] - other.m_data[0], m_data[1] - other.m_data[1], m_data[2] - other.m_data[2]);
}

// Cross Product
Vec3D Vec3D::operator*(const Vec3D& other) const {
    return Vec3D(
        m_data[1] * other.m_data[2] - m_data[2] * other.m_data[1],
        m_data[2] * other.m_data[0] - m_data[0] * other.m_data[2],
        m_data[0] * other.m_data[1] - m_data[1] * other.m_data[0]
    );
}

// Scalar Multiplication
Vec3D Vec3D::operator*(double scalar) const {
    return Vec3D(m_data[0] * scalar, m_data[1] * scalar, m_data[2] * scalar);
}

// Assignment operator
Vec3D& Vec3D::operator=(const Vec3D& other) {
    if (this != &other) { // Prevent self-assignment
        // Copy each element individually
        m_data[0] = other.m_data[0];
        m_data[1] = other.m_data[1];
        m_data[2] = other.m_data[2];
    }
    return *this; // Return *this to allow chained assignments
}
// Normalize Vector
void Vec3D::normalize() {
    double n = std::sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
    if (n > 0) {
        double inv_n = 1.0 / n; // Pre-calculate the reciprocal
        m_data[0] *= inv_n;
        m_data[1] *= inv_n;
        m_data[2] *= inv_n;
    }
}

// Vector Norm
double Vec3D::norm() const {
    return std::sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
}

// Dot Product
double Vec3D::dot(const Vec3D& v1, const Vec3D& v2) {
    return v1.m_data[0] * v2.m_data[0] + v1.m_data[1] * v2.m_data[1] + v1.m_data[2] * v2.m_data[2];
}

// Check if values are invalid
bool Vec3D::isbad() const {
    return !std::isfinite(m_data[0]) || !std::isfinite(m_data[1]) || !std::isfinite(m_data[2]);
}

// Check if values are valid
bool Vec3D::isgood() const {
    return std::isfinite(m_data[0]) && std::isfinite(m_data[1]) && std::isfinite(m_data[2]);
}

// Print Vector
void Vec3D::print() const {
    std::cout << m_data[0] << " " << m_data[1] << " " << m_data[2] << "\n";
}

// Stream Output
std::ostream& operator<<(std::ostream& os, const Vec3D& vec) {
    os << vec.m_data[0] << " " << vec.m_data[1] << " " << vec.m_data[2];
    return os;
}
