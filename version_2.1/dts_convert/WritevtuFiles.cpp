


#include "WritevtuFiles.h"


WritevtuFiles::WritevtuFiles()
{


}
WritevtuFiles::~WritevtuFiles()
{
    
}
void WritevtuFiles::Write(MeshBluePrint blueprint, std::string Filename)
{
    
    Vec3D Box=blueprint.simbox;
    std::vector<Vertex_Map> bvertex = blueprint.bvertex;
    std::vector<Triangle_Map> btriangle = blueprint.btriangle;
    std::vector<Inclusion_Map> binclusion = blueprint.binclusion;
    int numv=bvertex.size();
    int numtri=btriangle.size();
    int numtrirep=0;
    
    
    for (std::vector<Triangle_Map>::iterator it = btriangle.begin() ; it != btriangle.end(); ++it)
    {
        if(CrossPBC(Box, bvertex.at(it->v1), bvertex.at(it->v2), bvertex.at(it->v3))==false)
            numtrirep++;
    }
    std::ofstream Output;
    Output.open(Filename.c_str());
    

    
    
    
    
    Output<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<"\n";
    Output<<"  <UnstructuredGrid>"<<"\n";
    Output<<"    <Piece NumberOfPoints=\""<<numv<<"\" NumberOfCells=\""<<numtrirep<<"\">"<<"\n";
    Output<<"      <PointData Scalars=\"scalars\">"<<"\n";
    Output<<"        <DataArray type=\"Float32\" Name=\"inc\" Format=\"ascii\">"<<"\n";

    
    std::vector<int> vinc;
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
        vinc.push_back(0);
    
    for (std::vector<Inclusion_Map>::iterator it = binclusion.begin() ; it != binclusion.end(); ++it)
    {
        vinc.at(it->vid) = 1;
    }
    for (int i=0;i<numv;i++)
        Output<<"          "<<vinc.at(i)<<"\n";
        Output<<"        </DataArray>"<<"\n";
    Output<<"      </PointData>"<<"\n";
    Output<<"      <Points>"<<"\n";
    Output<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">"<<"\n";
    
    
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
        Output<<"          "<<it->x<<" "<<it->y<<" "<<it->z<<" "<<"\n";

    Output<<"        </DataArray>"<<"\n";
    Output<<"      </Points>"<<"\n";
    Output<<"      <Cells>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">"<<"\n";
    
    
    for (std::vector<Triangle_Map>::iterator it = btriangle.begin() ; it != btriangle.end(); ++it)
    {
        if(CrossPBC(Box, bvertex.at(it->v1), bvertex.at(it->v2), bvertex.at(it->v3))==false)
            Output<<"           "<<it->v1<<" "<<it->v2<<" "<<it->v3<<" "<<"\n";
    }
    Output<<"        </DataArray>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">"<<"\n";
    int ofset=3;
    Output<<"          ";
    for (int i=0;i<numtrirep;i++)
        Output<<ofset+3*i<<" ";
    
    Output<<"\n";
    Output<<"        </DataArray>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">"<<"\n";
    Output<<"          ";
    for (int i=0;i<numtrirep;i++)
        Output<<5<<" ";

    Output<<"\n";
    Output<<"        </DataArray>"<<"\n";
    Output<<"      </Cells>"<<"\n";
    Output<<"    </Piece>"<<"\n";
    Output<<"  </UnstructuredGrid>"<<"\n";
    Output<<"</VTKFile> "<<"\n";

}
bool WritevtuFiles::CrossPBC(Vec3D Box, Vertex_Map v1, Vertex_Map v2, Vertex_Map v3)
{
    bool is = false;
    
    Vec3D x1(v1.x,v1.y,v1.z);
    Vec3D x2(v2.x,v2.y,v2.z);
    Vec3D x3(v3.x,v3.y,v3.z);
    
    Vec3D X1 = x1-x2;
    Vec3D X2 = x1-x3;
    Vec3D X3 = x3-x2;
    
    for (int i=0;i<3;i++)
    {
        if(fabs(X1(i))>Box(i)/2.0)
            is=true;
        if(fabs(X2(i))>Box(i)/2.0)
            is=true;
        if(fabs(X3(i))>Box(i)/2.0)
            is=true;
    }
    return is;
}
