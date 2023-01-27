#if !defined(AFX_LMatrix_H_7F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_)
#define AFX_LMatrix_H_7F4421B8_D12D_11D3_CF24_124095086555__INCLUDED_

#include <vector>

/*

3 5 2
2 1 4    M=2 N=3



*/
class LMatrix
{


public:
    
	LMatrix(int M, int N);
    LMatrix(int M, int N, char );
	 ~LMatrix();
    
public:
    inline const int         Rsize()        const {return m_M;}
    inline const int         Csize()        const {return m_N;}  
    inline const std::vector <std::vector <double> >         MAT()        const {return m_Mat;} 

private:
int m_N;
int m_M;  
std::vector <std::vector <double> > m_Mat;

double m_T;
private:


public:
double norm();
double at( int n, int m);
void put( int n, int m, double);
void print();

LMatrix Inverse (LMatrix M);
LMatrix TM (LMatrix M);
double& operator()(const int n, const int m);
LMatrix operator + (LMatrix);
LMatrix operator - (LMatrix);
LMatrix operator * (LMatrix);
void operator = (LMatrix);

};

#endif
