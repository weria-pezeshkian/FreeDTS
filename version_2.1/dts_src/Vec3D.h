#ifndef VEC3D_H
#define VEC3D_H
/**
 * @file Vec3D.h
 * @brief A class representing a 3D vector with various mathematical operations.
 *
 * The Vec3D class provides a representation of 3-dimensional vectors and implements
 * essential vector operations such as addition, subtraction, scalar multiplication,
 * normalization, and dot product. This class is particularly useful in fields such
 * as physics, computer graphics, and computational geometry, where 3D vector
 * calculations are frequent.
 *
 * Features:
 * - Constructor for initializing a vector with specified coordinates or defaulting
 *   to the origin (0, 0, 0).
 * - Overloaded operators for intuitive vector arithmetic, including vector addition,
 *   subtraction, and cross product.
 * - Methods for computing the vector's norm (magnitude) and normalization.
 * - Static method to calculate the dot product between two Vec3D objects.
 * - Methods to check the validity of vector components, ensuring they are finite
 *   values.
 * - Support for outputting the vector to standard output via an overloaded
 *   `operator<<`.
 *
 * Usage Example:
 * - Create a vector: `Vec3D vec(1.0, 2.0, 3.0);`
 * - Perform vector operations: `Vec3D result = vec1 + vec2;`
 * - Normalize a vector: `vec.normalize();`
 * - Print the vector: `std::cout << vec;`
 *
 * @note This class uses an array for internal storage of vector components, which
 *       allows for fast indexing and efficient memory usage.
 */

#include <iostream>
#include <cmath>

class Vec3D {


public:
    Vec3D(double x = 0.0, double y = 0.0, double z = 0.0);
    
    inline double& operator()(int n) { return m_data[n]; }
    inline const double& operator()(int n) const { return m_data[n]; }

    Vec3D operator+(const Vec3D&) const;
    Vec3D operator-(const Vec3D&) const;
    Vec3D operator*(const Vec3D&) const;  // Cross Product
    Vec3D operator*(double) const;
    Vec3D& operator=(const Vec3D& other);

    void normalize();
    double norm() const;
    
    static double dot(const Vec3D&, const Vec3D&);

    bool isbad() const;
    bool isgood() const;
    void print() const;

    friend std::ostream& operator<<(std::ostream&, const Vec3D&);
    
private:
    double m_data[3];  // Array-based storage for fast indexing
};

#endif
