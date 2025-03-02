#include "Tensor2.h"

// Constructor
Tensor2::Tensor2(const Vec3D& x, const Vec3D& y, const Vec3D& z) {
    m_matrix[0][0] = x(0); m_matrix[0][1] = x(1); m_matrix[0][2] = x(2);
    m_matrix[1][0] = y(0); m_matrix[1][1] = y(1); m_matrix[1][2] = y(2);
    m_matrix[2][0] = z(0); m_matrix[2][1] = z(1); m_matrix[2][2] = z(2);
}

// Default Constructor
Tensor2::Tensor2() {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m_matrix[i][j] = 0.0; // Initialize all elements to zero
        }
    }
}

// Special Matrix Constructor
Tensor2::Tensor2(char t) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            m_matrix[i][j] = 0.0; // Initialize all to zero
        }
    }
    if (t == 'I') {
        m_matrix[0][0] = 1.0; m_matrix[1][1] = 1.0; m_matrix[2][2] = 1.0; // Identity matrix
    }
}

// Destructor
Tensor2::~Tensor2() {}

// Access element
double Tensor2::at(int n, int m) const {
    return m_matrix[n][m]; // Direct access without out of range checks
}

// Put element
void Tensor2::put(int n, int m, double s) {
    m_matrix[n][m] = s; // Direct access without out of range checks
}

// Overloaded operator for element access
double& Tensor2::operator()(int n, int m) {
    return m_matrix[n][m]; // Direct access
}

// Matrix-vector multiplication
Vec3D Tensor2::operator*(const Vec3D& A) const {
    return Vec3D(
        m_matrix[0][0] * A(0) + m_matrix[0][1] * A(1) + m_matrix[0][2] * A(2),
        m_matrix[1][0] * A(0) + m_matrix[1][1] * A(1) + m_matrix[1][2] * A(2),
        m_matrix[2][0] * A(0) + m_matrix[2][1] * A(1) + m_matrix[2][2] * A(2)
    );
}

// Matrix addition
Tensor2 Tensor2::operator+(const Tensor2& M) const {
    Tensor2 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result(i, j) = this->m_matrix[i][j] + M.m_matrix[i][j];
        }
    }
    return result;
}

// Matrix subtraction
Tensor2 Tensor2::operator-(const Tensor2& M) const {
    Tensor2 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result(i, j) = this->m_matrix[i][j] - M.m_matrix[i][j];
        }
    }
    return result;
}

// Matrix multiplication
Tensor2 Tensor2::operator*(const Tensor2& M) const {
    Tensor2 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result(i, j) = m_matrix[i][0] * M.m_matrix[0][j] +
                           m_matrix[i][1] * M.m_matrix[1][j] +
                           m_matrix[i][2] * M.m_matrix[2][j];
        }
    }
    return result;
}

// Matrix scalar multiplication
Tensor2 Tensor2::operator*(double x) const {
    Tensor2 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result(i, j) = this->m_matrix[i][j] * x;
        }
    }
    return result;
}

// Transpose the matrix
Tensor2 Tensor2::Transpose() const {
    Tensor2 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result(j, i) = this->m_matrix[i][j];
        }
    }
    return result;
}

// Create a tensor from a vector
Tensor2 Tensor2::makeTen(const Vec3D& X) {
    Tensor2 A;

     // Precompute the components of X
     const double x0 = X(0);
     const double x1 = X(1);
     const double x2 = X(2);

     // Manually unroll the loop for better performance
     A(0, 0) = x0 * x0; A(0, 1) = x0 * x1; A(0, 2) = x0 * x2;
     A(1, 0) = x1 * x0; A(1, 1) = x1 * x1; A(1, 2) = x1 * x2;
     A(2, 0) = x2 * x0; A(2, 1) = x2 * x1; A(2, 2) = x2 * x2;

     return A;
}
