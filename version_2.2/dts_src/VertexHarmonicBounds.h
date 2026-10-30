#ifndef VertexHarmonicBounds_H
#define VertexHarmonicBounds_H

#include "SimDef.h"
#include "Vec3D.h"
#include "bond.h"


class MESH;
class VertexHarmonicBounds {

public:

    VertexHarmonicBounds();
    ~VertexHarmonicBounds();
    inline std::vector<bond *> GetBonds()   {return m_pVertexBond;}
    double GetBondEnergyOfVertex();
    void AddBondToList(bond* b);

protected:
    std::vector<bond *> m_pVertexBond;


};

#endif // VertexHarmonicBounds_H
