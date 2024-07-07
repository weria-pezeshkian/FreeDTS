#if !defined(AFX_Vec3D_H_8F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_)
#define AFX_Vec3D_H_8F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_

#include <iostream>
#include <cmath>
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An inhouse made  3d vector object.
 *******************/
//class Vec3D;
//std::ostream& operator<<(std::ostream& os, const Vec3D& vec);

class Vec3D {
public:
    Vec3D(double x = 0.0, double y = 0.0, double z = 0.0);
    double& operator()(const int n);
    Vec3D operator+(const Vec3D& other) const;
    Vec3D operator-(const Vec3D& other) const;
    Vec3D operator*(const Vec3D& other) const;
    Vec3D operator*(double scalar) const;
    void operator=(const Vec3D& other);
    double norm() const;
    static double dot(const Vec3D& v1, const Vec3D& v2);
    bool isbad() const; // Function to check if any element is non-finite
    bool isgood() const;
    void print();
    void normalize();
    friend std::ostream& operator<<(std::ostream& os, const Vec3D& vec);


private:
    double m_X;
    double m_Y;
    double m_Z;
};

#endif
