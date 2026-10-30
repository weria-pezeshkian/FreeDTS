#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
#include <cmath>

class Vec3D {
public:
    // Make ALL constructors constexpr
    constexpr Vec3D(double x = 0.0, double y = 0.0, double z = 0.0) noexcept
        : m_data{x, y, z} {}
    
    constexpr Vec3D(const Vec3D& other) noexcept = default;
    constexpr Vec3D(Vec3D&& other) noexcept = default;
    
    constexpr double& operator()(int n) noexcept { return m_data[n]; }
    constexpr const double& operator()(int n) const noexcept { return m_data[n]; }

    // These can't be constexpr because they return *this (but that's fine)
    Vec3D& operator=(const Vec3D& other) noexcept = default;
    Vec3D& operator=(Vec3D&& other) noexcept = default;

    // Arithmetic operators - make constexpr
    constexpr Vec3D operator+(const Vec3D& other) const noexcept {
        return Vec3D(m_data[0] + other.m_data[0],
                     m_data[1] + other.m_data[1],
                     m_data[2] + other.m_data[2]);
    }
    
    constexpr Vec3D operator-(const Vec3D& other) const noexcept {
        return Vec3D(m_data[0] - other.m_data[0],
                     m_data[1] - other.m_data[1],
                     m_data[2] - other.m_data[2]);
    }
    
    // Cross product
    constexpr Vec3D operator*(const Vec3D& other) const noexcept {
        return Vec3D(
            m_data[1] * other.m_data[2] - m_data[2] * other.m_data[1],
            m_data[2] * other.m_data[0] - m_data[0] * other.m_data[2],
            m_data[0] * other.m_data[1] - m_data[1] * other.m_data[0]
        );
    }
    
    // Scalar multiplication
    constexpr Vec3D operator*(double scalar) const noexcept {
        return Vec3D(m_data[0] * scalar,
                     m_data[1] * scalar,
                     m_data[2] * scalar);
    }

    // These can't be constexpr (use sqrt, modify this)
    void normalize() noexcept {
        double n = std::sqrt(m_data[0] * m_data[0] +
                            m_data[1] * m_data[1] +
                            m_data[2] * m_data[2]);
        if (n > 0) {
            double inv_n = 1.0 / n;
            m_data[0] *= inv_n;
            m_data[1] *= inv_n;
            m_data[2] *= inv_n;
        }
    }
    
    double norm() const noexcept {
        return std::sqrt(m_data[0] * m_data[0] +
                         m_data[1] * m_data[1] +
                         m_data[2] * m_data[2]);
    }
    
    // constexpr version without sqrt
    constexpr double norm_squared() const noexcept {
        return m_data[0] * m_data[0] +
               m_data[1] * m_data[1] +
               m_data[2] * m_data[2];
    }
    
    // Static dot product - constexpr
    static constexpr double dot(const Vec3D& v1, const Vec3D& v2) noexcept {
        return v1.m_data[0] * v2.m_data[0] +
               v1.m_data[1] * v2.m_data[1] +
               v1.m_data[2] * v2.m_data[2];
    }

    constexpr bool isbad() const noexcept {
        // Note: std::isfinite is NOT constexpr in C++17
        // So this can't be constexpr with the check
        return !std::isfinite(m_data[0]) ||
               !std::isfinite(m_data[1]) ||
               !std::isfinite(m_data[2]);
    }
    
    constexpr bool isgood() const noexcept {
        return !isbad();
    }
    
    void print() const {
        std::cout << m_data[0] << ' ' << m_data[1] << ' ' << m_data[2] << '\n';
    }

    friend std::ostream& operator<<(std::ostream& os, const Vec3D& vec) {
        os << vec.m_data[0] << ' ' << vec.m_data[1] << ' ' << vec.m_data[2];
        return os;
    }
    
private:
    double m_data[3];
};

#endif
