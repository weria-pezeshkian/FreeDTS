#if !defined(AFX_AnisotropicStretchingByNematicField_H_334B21B8_INCLUDED_)
#define AFX_AnisotropicStretchingByNematicField_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractLocalStretching.h"
#include "AbstractForceonVerticesfromInclusions.h"
class AnisotropicStretchingByNematicField : public AbstractLocalStretching{
public:
    AnisotropicStretchingByNematicField(double K_p, double K_n);
    ~AnisotropicStretchingByNematicField();
    double Energy(vertex *p);
    std::string CurrentState();
    void Initialize();
    inline  std::string GetDerivedDefaultReadName()  {return "AnisotropicStretchingByNematicField";}
    inline static std::string GetDefaultReadName() {return "AnisotropicStretchingByNematicField";}

private:


    double m_Kp;
    double m_Kn;

    

};


#endif
