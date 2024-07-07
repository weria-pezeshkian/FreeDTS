#if !defined(TEM_Voxel_H_INCLUDED_)
#define TEM_Voxel_H_INCLUDED_

#include "SimDef.h"
#include "Voxelization.h"
/*
    Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
    Updated: 2024

    Description:
    The Voxel class represents a single voxel within a 3D grid used for voxelizing a box-shaped region. It inherits from the Voxelization class.

    Voxel Characteristics:
    - Each voxel contains a list of objects (vertices/points) that fall within its bounds.
    - The voxel grid is dynamic, adjusting its size and contents during simulation time.

    Dependencies:
    - SimDef.h: Contains definitions and constants used in the simulation.
    - Voxelization.h: Declaration of the Voxelization class, which handles the voxelization process.

    Usage Example:
    // Include necessary headers
    #include "Voxel.h"
    #include "SomeObjectType.h"

    // Create an instance of Voxel with appropriate parameters
    Voxel<SomeObjectType> voxel();

    // Access voxel properties and functionalities
    voxel.AddtoContentList();

    Class Template:
    - Voxel is templated on Type, allowing it to work with different types of objects.

    Important Notes:
    - The Voxel class assumes ownership of the objects added to its content list. Ensure proper memory management when removing objects.
*/
template<typename Type>
class Voxel  {
public:
    Voxel(int id, int n, int m, int k, int &nx, int &ny, int &nz)
        : m_ID(id), m_X_index(n), m_Y_index(m), m_Z_index(k),
          m_Nx(nx), m_Ny(ny), m_Nz(nz) {
        m_NeighbouringVoxel[1][1][1] = this;
    }
    ~Voxel(){
        m_List.clear();
    }// ~Voxel()

    inline const int          GetID()       const        { return m_ID; }      // voxel id
    inline int                GetXIndex()                { return m_X_index; } // x index,  m_X_index*m_Vox_Lx 
    inline int                GetYIndex()                { return m_Y_index; } //
    inline int                GetZIndex()                { return m_Z_index; } //
    inline double GetXSideVoxel(double Lx)     const { return Lx/double(m_Nx); }
    inline double GetYSideVoxel(double Ly)     const { return Ly/double(m_Ny); }
    inline double GetZSideVoxel(double Lz)     const { return Lz/double(m_Nz); }
    inline const int GetXNoVoxel()                   const { return m_Nx; }
    inline const int GetYNoVoxel()                   const { return m_Ny; }
    inline const int GetZNoVoxel()                   const { return m_Nz; }
    inline std::vector<Type*> GetContentObjects()        { return m_List; }    //all the objects, e.g., vertices in the voxel

public:
    void AddtoContentList(Type * p_obj){ // adding a new vertex/point, when the vertex enters the domain
        m_List.push_back(p_obj);
        return;
    }
    void ClearList(){ // adding a new vertex/point, when the vertex enters the domain
        m_List.clear();
        return;
    }
    bool SetANeighbourCell(int i, int j, int k, Voxel *vox){// setting the nb voxels 26+1 voxels
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            std::cout << " ---> error: Invalid indices. Indices must be within range [-1, 1]." << std::endl;
            return false;
        }

        m_NeighbouringVoxel[i + 1][j + 1][k + 1] = vox;
        return true;
    }
    void RemoveObjectFromContentList(Type * p_obj){// removing a  vertex, when the vertex exist the domain
        m_List.erase(std::remove(m_List.begin(), m_List.end(), p_obj), m_List.end());
        return;
    }
    Voxel * GetANeighbourCell(int i, int j, int k){// accessing a specific nighbouring cell
        if(i==0 && j==0 && k==0){
            return this;
        }
        else if(i>1 || i<-1 || j>1 || j<-1 || k>1 || k<-1  )
            std::cout<<" ---> error: such indices are not permitted "<<std::endl;
        
        return m_NeighbouringVoxel[i+1][j+1][k+1];
    }
    static int Convert2LocalVoxelIndex(int new_index, int old_index, int No_index){

        if(new_index == old_index){
            return 0;
        }
        else if( (old_index + 1 + No_index) % No_index == new_index){
            return 1;
        }
        else if( (old_index - 1 + No_index) % No_index == new_index){
            return -1;
        }
        else {
            std::cout << " ---> error. the object has moved  out of neighbouring voxel " << std::endl;
            std::cout <<"new_id: " <<new_index <<" old_id: "<< old_index<<" no-voxel: "<<No_index<< std::endl;
        }
        return 0;
    }
private:
  int m_ID;
  std::vector <Type *> m_List;              // List of the vertices/points that the cell contains/ Content
    int m_X_index;                          // indices of the cell withing the box;
    int m_Y_index;
    int m_Z_index;
    Voxel *m_NeighbouringVoxel[3][3][3];
    int &m_Nx; // Number of the voxels in the x direction
    int &m_Ny; // Number of the voxels in the y direction
    int &m_Nz; // Number of the voxels in the z direction


};

#endif
