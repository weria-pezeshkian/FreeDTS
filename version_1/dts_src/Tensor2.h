#if !defined(AFX_Tensor2_H_3O4421B8_D12D_11D3_CF24_124095086555__INCLUDED_)
#define AFX_Tensor2_H_3O4421B8_D12D_11D3_CF24_124095086555__INCLUDED_
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An inhouse  3x3 matrix object.
 This can be made more effiently!!!! checking later
 *******************/
#include <vector>
#include "Vec3D.h"
class Tensor2
{
public:
	Tensor2(Vec3D v1,Vec3D v2,Vec3D v3);
    Tensor2();
    Tensor2(char);
	 ~Tensor2();
public:
    inline Vec3D         GetV1()        const {return m_V1;}
    inline Vec3D         GetV2()        const {return m_V2;}
    inline Vec3D         GetV3()        const {return m_V3;}
    
private:
    Vec3D m_V1;
    Vec3D m_V2;
    Vec3D m_V3;

public:
double at( int n,int m);
void put( int n, int m, double);
Vec3D MULT(Tensor2, Vec3D);


double& operator()(const int n,const int);
Tensor2 operator + (Tensor2);
Tensor2 operator - (Tensor2);
    
Tensor2 operator * (Tensor2);
Tensor2 operator * (double);
Vec3D  operator * (Vec3D);

void operator = (Tensor2);
    Tensor2 makeTen(Vec3D);
    Tensor2 Transpose(Tensor2);

};

#endif
