#ifndef VOXELIZATION_H_INCLUDED
#define VOXELIZATION_H_INCLUDED
#ifdef _OPENMP
# include <omp.h>
#endif
#include "SimDef.h"
#include "Vec3D.h"
/*
    Author: Weria Pezeshkian (2024)

    Description:
    The Voxelization class divides a box-shaped region into a grid of voxels and assigns objects (such as vertices or points) to these voxels based on their positions. It provides functionality for voxelizing objects, updating voxel sizes, and reallocating objects to voxels when necessary.

    Voxelization Process:
    - The box region is divided into voxels based on a specified voxel size.
    - Objects are assigned to voxels based on their positions within the box.
    - If the box dimensions change, so the size of a voxel is smaller then min distance,
      the voxelization can be updated by rerunning the Voxelize function.

    Dependencies:
    - SimDef.h: Contains definitions and constants used in the simulation.
    - Vec3D.h: Provides a class for representing 3D vectors.

    Usage Example:
    ```
    // Include necessary headers
    #include "Voxelization.h"
    #include "SomeObjectType.h"

    // Create an instance of Voxelization with appropriate parameters
    Voxelization<SomeObjectType> voxelization(...);

    // Perform voxelization by calling the Voxelize function
    voxelization.Voxelize(...);

    // Use voxelization results as needed

    Class Template:
    - Voxelization is templated on Type, allowing it to work with different types of objects.
 
    functional members:
      GetTSideVoxel()    T=X,Y,Z (3 functions)
      GetTVoxelNumber()  T=X,Y,Z (3 functions)
      UpdateVoxelSize()
      Voxelize()
      ReassignMembersToVoxels()
*/
template<typename Type>
class Voxel;       // Forward declaration of the Voxel class
template<typename Type>
class Voxelization
{
public:
Voxelization(){ // Use Type here
    m_Lx = 1;
    m_Ly = 1;
    m_Lz = 1;
}
Voxelization(Vec3D *pBox) { // Use Type here
        // Initialize the box object
        m_pBox = pBox;
        m_Lx = 1;
        m_Ly = 1;
        m_Lz = 1;

}
~Voxelization() {// to clear the memory allocated for this
        if (m_AllVoxel != nullptr) {
            for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                delete m_AllVoxel[i][j][k];
            }
                delete[] m_AllVoxel[i][j];
            }
                delete[] m_AllVoxel[i];
            }
            delete[] m_AllVoxel;
            m_AllVoxel = nullptr;
        }
}
public:
inline double GetXSideVoxel() const { return (*m_pBox)(0)/double(m_Nx); }
inline double GetYSideVoxel() const { return (*m_pBox)(1)/double(m_Ny); }
inline double GetZSideVoxel() const { return (*m_pBox)(2)/double(m_Nz); }
inline int GetXVoxelNumber() const { return m_Nx; }
inline int GetYVoxelNumber() const { return m_Ny; }
inline int GetZVoxelNumber() const { return m_Nz; }
inline  Voxel<Type> **** GetAllVoxel() const { return m_AllVoxel;}



void UpdateVoxelSize(double lx, double ly, double lz) {
    m_Lx = lx;
    m_Ly = ly;
    m_Lz = lz;
    return;
}
void SetBox(Vec3D *pBox) {
    m_pBox = pBox;
    return;
}
//----> Important and detailed function
bool Voxelize(std::vector<Type *> all_pObjects) {
#if DEBUG_MODE == Enabled
    std::cout<<" We are in the Voxelize function \n";
#endif
    
    //---> make m_AllVoxel empty
            if (m_AllVoxel != nullptr) {
                for (int i = 0; i < m_Nx; ++i) {
                for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    delete m_AllVoxel[i][j][k];
                }
                    delete[] m_AllVoxel[i][j];
                }
                    delete[] m_AllVoxel[i];
                }
                delete[] m_AllVoxel;
                m_AllVoxel = nullptr;
            }
#if DEBUG_MODE == Enabled
    std::cout<<" all voxels has been emptied \n";
#endif
        //--->find the appropriate voxel size and number of voxels
        m_Nx = int((*m_pBox)(0)/m_Lx);
        m_Ny = int((*m_pBox)(1)/m_Ly);
        m_Nz = int((*m_pBox)(2)/m_Lz);
#if DEBUG_MODE == Enabled
    std::cout<<m_Nx<<"  "<<m_Ny<<"  "<<m_Nz<<" number of voxels\n";
#endif
        double Lx = (*m_pBox)(0)/double(m_Nx);
        double Ly = (*m_pBox)(1)/double(m_Ny);
        double Lz = (*m_pBox)(2)/double(m_Nz);
#if DEBUG_MODE == Enabled
    std::cout<<Lx<<"  "<<Ly<<"  "<<Lz<<" size of voxels\n";
#endif
        //----> making m_AllVoxel an m_AllVoxel[m_Nx][m_Ny][m_Nz] and create all the cells
        int voxel_id = 0;
        m_AllVoxel = new Voxel<Type>***[m_Nx];
        for (int i = 0; i < m_Nx; ++i) {
            m_AllVoxel[i] = new Voxel<Type>**[m_Ny];
            for (int j = 0; j < m_Ny; ++j) {
                m_AllVoxel[i][j] = new Voxel<Type>*[m_Nz];
                for (int k = 0; k < m_Nz; ++k) {
                    m_AllVoxel[i][j][k] = new Voxel<Type>(voxel_id,i,j,k, m_Nx, m_Ny, m_Nz);
                    voxel_id++;
                }
            }
        }
#if DEBUG_MODE == Enabled
    std::cout<<" we created new voxels \n";
#endif
        //--> Update voxel neighbours
        //----> making m_AllVoxel an m_AllVoxel[m_Nx][m_Ny][m_Nz] and create all the cells
        for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    for (int n = -1; n <=1; ++n)
                    for (int m = -1; m <= 1; ++m)
                    for (int p = -1; p <= 1; ++p){
                        (m_AllVoxel[i][j][k])->SetANeighbourCell(n,m,p, m_AllVoxel[(i+n+m_Nx)%m_Nx][(j+m+m_Ny)%m_Ny][(k+p+m_Nz)%m_Nz]);
                    }
                }
            }
        }
#if DEBUG_MODE == Enabled
    std::cout<<" voxel nieghbours has been set \n";
#endif
        //--->allocate the objects to the voxels
        for (typename std::vector<Type *>::iterator it = all_pObjects.begin(); it != all_pObjects.end(); ++it) {
            double x = (*it)->GetXPos(); // Retrieve X position of each object
            double y = (*it)->GetYPos(); // Retrieve Y position of each object
            double z = (*it)->GetZPos(); // Retrieve Z position of each object
            // Process or allocate objects to voxels based on their XYZ position
            int nx = (*it)->GetXPos()/Lx;
            int ny = (*it)->GetYPos()/Ly;
            int nz = (*it)->GetZPos()/Lz;
            Voxel<Type> * pVoxel = m_AllVoxel[nx][ny][nz];
            (m_AllVoxel[nx][ny][nz])->AddtoContentList(*it);
            (*it)->UpdateVoxel(pVoxel);
        }
#if DEBUG_MODE == Enabled
    std::cout<<" vertices are added to voxels \n";
#endif
        return true; // Placeholder
}
bool VoxelizeOpenMP(std::vector<Type *> all_pObjects) {
#ifdef _OPENMP
        #if DEBUG_MODE == Enabled
        std::cout << " We are in the Voxelize function using openMP\n";
        #endif

        //---> Make m_AllVoxel empty
        if (m_AllVoxel != nullptr) {
            for (int i = 0; i < m_Nx; ++i) {
                for (int j = 0; j < m_Ny; ++j) {
                    for (int k = 0; k < m_Nz; ++k) {
                        delete m_AllVoxel[i][j][k];
                    }
                    delete[] m_AllVoxel[i][j];
                }
                delete[] m_AllVoxel[i];
            }
            delete[] m_AllVoxel;
            m_AllVoxel = nullptr;
        }
        #if DEBUG_MODE == Enabled
        std::cout << " All voxels have been emptied \n";
        #endif

        //---> Find the appropriate voxel size and number of voxels
        m_Nx = int((*m_pBox)(0) / m_Lx);
        m_Ny = int((*m_pBox)(1) / m_Ly);
        m_Nz = int((*m_pBox)(2) / m_Lz);
        double Lx = (*m_pBox)(0) / double(m_Nx);
        double Ly = (*m_pBox)(1) / double(m_Ny);
        double Lz = (*m_pBox)(2) / double(m_Nz);

        #if DEBUG_MODE == Enabled
        std::cout << m_Nx << "  " << m_Ny << "  " << m_Nz << " number of voxels\n";
        std::cout << Lx << "  " << Ly << "  " << Lz << " size of voxels\n";
        #endif

        //---> Create voxels
        int voxel_id = 0;
        m_AllVoxel = new Voxel<Type>***[m_Nx];

        // Parallelize voxel creation
    // First allocate memory sequentially to avoid race conditions
    for (int i = 0; i < m_Nx; ++i) {
        m_AllVoxel[i] = new Voxel<Type>**[m_Ny];
        for (int j = 0; j < m_Ny; ++j) {
            m_AllVoxel[i][j] = new Voxel<Type>*[m_Nz];
        }
    }
    // Now parallelize the voxel creation safely
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                int voxel_id = i * m_Ny * m_Nz + j * m_Nz + k;
                m_AllVoxel[i][j][k] = new Voxel<Type>(voxel_id, i, j, k, m_Nx, m_Ny, m_Nz);
            }
        }
    }
        #if DEBUG_MODE == Enabled
        std::cout << " We created new voxels \n";
        #endif

        //---> Update voxel neighbors
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    for (int n = -1; n <= 1; ++n) {
                        for (int m = -1; m <= 1; ++m) {
                            for (int p = -1; p <= 1; ++p) {
                                m_AllVoxel[i][j][k]->SetANeighbourCell(
                                    n, m, p,
                                    m_AllVoxel[(i + n + m_Nx) % m_Nx]
                                              [(j + m + m_Ny) % m_Ny]
                                              [(k + p + m_Nz) % m_Nz]
                                );
                            }
                        }
                    }
                }
            }
        }

        #if DEBUG_MODE == Enabled
        std::cout << " Voxel neighbors have been set \n";
        #endif

        //---> Allocate objects to the voxels
        #pragma omp parallel for
        for (typename std::vector<Type *>::iterator it = all_pObjects.begin(); it != all_pObjects.end(); ++it) {
            double x = (*it)->GetXPos(); // Retrieve X position of each object
            double y = (*it)->GetYPos(); // Retrieve Y position
            double z = (*it)->GetZPos(); // Retrieve Z position

            int nx = x / Lx;
            int ny = y / Ly;
            int nz = z / Lz;

            Voxel<Type> *pVoxel = m_AllVoxel[nx][ny][nz];

            // Critical section to avoid race conditions when adding objects to voxel content
            #pragma omp critical
            {
                pVoxel->AddtoContentList(*it);
                (*it)->UpdateVoxel(pVoxel);
            }
        }

        #if DEBUG_MODE == Enabled
        std::cout << " Objects are added to voxels \n";
        #endif

        return true; // Placeholder
#else
    std::cerr << "---> Error: OpenMP is not available, but the program requires it. Please recompile with OpenMP support.\n";
    exit(EXIT_FAILURE); // Exit with a non-zero value to indicate failure
    
    return false;
#endif
}
//----> Important and detailed function
bool ReassignMembersToVoxels(std::vector<Type *> all_pObjects) {
            
            double Lx = (*m_pBox)(0)/double(m_Nx);
            double Ly = (*m_pBox)(1)/double(m_Ny);
            double Lz = (*m_pBox)(2)/double(m_Nz);
    //---> remove all the the members from the voxels
            if (m_AllVoxel != nullptr) {
                for (int i = 0; i < m_Nx; ++i) {
                for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    m_AllVoxel[i][j][k]->ClearList();
                }}}
            } //if (m_AllVoxel != nullptr) {
            //--->allocate the objects to the voxels
            for (typename std::vector<Type *>::iterator it = all_pObjects.begin(); it != all_pObjects.end(); ++it) {
                double x = (*it)->GetXPos(); // Retrieve X position of each object
                double y = (*it)->GetYPos(); // Retrieve Y position of each object
                double z = (*it)->GetZPos(); // Retrieve Z position of each object
                // Process or allocate objects to voxels based on their Z position
                int nx = (*it)->GetXPos()/Lx;
                int ny = (*it)->GetYPos()/Ly;
                int nz = (*it)->GetZPos()/Lz;
                Voxel<Type> * pVoxel = m_AllVoxel[nx][ny][nz];
                (m_AllVoxel[nx][ny][nz])->AddtoContentList(*it);
                (*it)->UpdateVoxel(pVoxel);
            }

 return true; // Placeholder
} // ReassignMembersToVoxels(std::vector<Type *> all_pObjects) {
private:
    int m_Nx; // Number of the voxels in the x direction
    int m_Ny; // Number of the voxels in the y direction
    int m_Nz; // Number of the voxels in the z direction
    double m_Lx; // voxel length in the x direction
    double m_Ly; // voxel length in the y direction
    double m_Lz; // voxel length in the z direction

    Vec3D *m_pBox;
 //   Voxel<Type> *m_AllVoxel[2][2][2];
    Voxel<Type> ****m_AllVoxel;

public:
    // Add other member functions here
};

#endif // VOXELIZATION_H_INCLUDED
