#ifndef TENSOR2_H
#define TENSOR2_H

#include "Vec3D.h"
#include <iostream>

class Tensor2 {
public:
    // Constructors - all constexpr
    constexpr Tensor2(const Vec3D& x, const Vec3D& y, const Vec3D& z) noexcept
        : m_matrix{
            x(0), x(1), x(2),
            y(0), y(1), y(2),
            z(0), z(1), z(2)
          } {}
    
    constexpr Tensor2() noexcept : m_matrix{} {}
    
    constexpr explicit Tensor2(char t) noexcept : m_matrix{} {
        if (t == 'I') {
            m_matrix[0] = 1.0;   // [0][0]
            m_matrix[4] = 1.0;   // [1][1]
            m_matrix[8] = 1.0;   // [2][2]
        }
    }
    
    ~Tensor2() = default;
    
    // Copy/Move operations
    constexpr Tensor2(const Tensor2&) noexcept = default;
    constexpr Tensor2(Tensor2&&) noexcept = default;
    constexpr Tensor2& operator=(const Tensor2&) noexcept = default;
    constexpr Tensor2& operator=(Tensor2&&) noexcept = default;
    
    // Access functions - using flat array for better performance
    constexpr double at(int n, int m) const noexcept {
        return m_matrix[n * 3 + m];
    }
    
    constexpr void put(int n, int m, double s) noexcept {
        m_matrix[n * 3 + m] = s;
    }
    
    constexpr double& operator()(int n, int m) noexcept {
        return m_matrix[n * 3 + m];
    }
    
    constexpr const double& operator()(int n, int m) const noexcept {
        return m_matrix[n * 3 + m];
    }
    
    // Matrix-vector multiplication
    constexpr Vec3D operator*(const Vec3D& A) const noexcept {
        return Vec3D(
            m_matrix[0] * A(0) + m_matrix[1] * A(1) + m_matrix[2] * A(2),
            m_matrix[3] * A(0) + m_matrix[4] * A(1) + m_matrix[5] * A(2),
            m_matrix[6] * A(0) + m_matrix[7] * A(1) + m_matrix[8] * A(2)
        );
    }
    
    // Matrix addition
    constexpr Tensor2 operator+(const Tensor2& M) const noexcept {
        Tensor2 result;
        for (int i = 0; i < 9; ++i) {
            result.m_matrix[i] = m_matrix[i] + M.m_matrix[i];
        }
        return result;
    }
    
    // Matrix subtraction
    constexpr Tensor2 operator-(const Tensor2& M) const noexcept {
        Tensor2 result;
        for (int i = 0; i < 9; ++i) {
            result.m_matrix[i] = m_matrix[i] - M.m_matrix[i];
        }
        return result;
    }
    
    // Matrix multiplication
    constexpr Tensor2 operator*(const Tensor2& M) const noexcept {
        Tensor2 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int idx = i * 3 + j;
                result.m_matrix[idx] =
                    m_matrix[i * 3 + 0] * M.m_matrix[0 * 3 + j] +
                    m_matrix[i * 3 + 1] * M.m_matrix[1 * 3 + j] +
                    m_matrix[i * 3 + 2] * M.m_matrix[2 * 3 + j];
            }
        }
        return result;
    }
    
    // Scalar multiplication
    constexpr Tensor2 operator*(double x) const noexcept {
        Tensor2 result;
        for (int i = 0; i < 9; ++i) {
            result.m_matrix[i] = m_matrix[i] * x;
        }
        return result;
    }
    
    // Transpose
    constexpr Tensor2 Transpose() const noexcept {
        Tensor2 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m_matrix[j * 3 + i] = m_matrix[i * 3 + j];
            }
        }
        return result;
    }
    
    // Create dyadic product
    static constexpr Tensor2 makeTen(const Vec3D& X) noexcept {
        double x0 = X(0), x1 = X(1), x2 = X(2);
        
        Tensor2 result;
        result.m_matrix[0] = x0 * x0;  // [0][0]
        result.m_matrix[1] = x0 * x1;  // [0][1]
        result.m_matrix[2] = x0 * x2;  // [0][2]
        result.m_matrix[3] = x1 * x0;  // [1][0]
        result.m_matrix[4] = x1 * x1;  // [1][1]
        result.m_matrix[5] = x1 * x2;  // [1][2]
        result.m_matrix[6] = x2 * x0;  // [2][0]
        result.m_matrix[7] = x2 * x1;  // [2][1]
        result.m_matrix[8] = x2 * x2;  // [2][2]
        return result;
    }
    
    // Utility functions
    constexpr double trace() const noexcept {
        return m_matrix[0] + m_matrix[4] + m_matrix[8];
    }
    
    constexpr double determinant() const noexcept {
        return m_matrix[0] * (m_matrix[4] * m_matrix[8] - m_matrix[5] * m_matrix[7])
             - m_matrix[1] * (m_matrix[3] * m_matrix[8] - m_matrix[5] * m_matrix[6])
             + m_matrix[2] * (m_matrix[3] * m_matrix[7] - m_matrix[4] * m_matrix[6]);
    }
    
    void print() const {
        for (int i = 0; i < 3; ++i) {
            std::cout << m_matrix[i * 3 + 0] << ' '
                      << m_matrix[i * 3 + 1] << ' '
                      << m_matrix[i * 3 + 2] << '\n';
        }
    }
    
    friend constexpr Tensor2 operator*(double x, const Tensor2& M) noexcept {
        return M * x;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Tensor2& T) {
        for (int i = 0; i < 3; ++i) {
            os << T.m_matrix[i * 3 + 0] << ' '
               << T.m_matrix[i * 3 + 1] << ' '
               << T.m_matrix[i * 3 + 2];
            if (i < 2) os << '\n';
        }
        return os;
    }
    
private:
    // Flat array for better cache locality and constexpr support
    double m_matrix[9]{};
};

#endif
