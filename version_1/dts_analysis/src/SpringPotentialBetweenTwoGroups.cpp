


#include "SpringPotentialBetweenTwoGroups.h"
#include "Nfunction.h"
SpringPotentialBetweenTwoGroups::SpringPotentialBetweenTwoGroups()
{
    m_State = false;
    m_K = 0;
    m_State = 0;
    m_L0 = 0;
    m_T0 = 0;

}
SpringPotentialBetweenTwoGroups::SpringPotentialBetweenTwoGroups(bool state)
{
    m_State = state;
    m_K = 0;
    m_State = 0;
    m_L0 = 0;
    m_T0 = 0;
}
SpringPotentialBetweenTwoGroups::SpringPotentialBetweenTwoGroups(bool state,  double K, double lt, int t0, std::string group1,std::string group2,double nx,double ny,double nz)
{
    m_K = K/2;
    m_State = state;
    m_Group1Name = group1;
    m_Group2Name = group2;
    m_Direction(0) = nx;
    m_Direction(1) = ny;
    m_Direction(2) = nz;
    m_LT = lt;
    m_T0 = t0;
    
    
    // We should prepare the system here
    //1) Make groups and find there com
    // create a function to update group com if a vertex moves, it should check if the vertex exist in the group and if yes, update the com.
    // We should also reject the move or accept it, need more functions
}

SpringPotentialBetweenTwoGroups::~SpringPotentialBetweenTwoGroups()
{
    
}
void SpringPotentialBetweenTwoGroups::CalculateEnergy(int step)
{
    double en=0;
    double f= 0;
    double l0 = 0;
    if(step>m_T0)
        l0= m_LT;
    else
        l0 = m_L0+(m_LT-m_L0)*step/m_T0;

    Vec3D *pBox = (m_pGroup1.at(0))->GetBox();
    Vec3D PX1 = m_Group2COG;
    Vec3D PX2 = m_Group2COG;
    Vec3D Dist;

    for (int i=0;i<3;i++)
    {
        double dist=PX1(i)-PX2(i);
        if(dist>(*pBox)(i)/2.0)
        {
            dist=(*pBox)(i)-dist;
        }
        else if(dist<0)
        {
            dist=-dist;
        }
        Dist(i)=dist;
    }
    double dist=m_Direction(0)*Dist(0)*Dist(0)+m_Direction(1)*Dist(1)*Dist(1)+m_Direction(2)*Dist(2)*Dist(2);
    dist=sqrt(dist);
    
    double m_Energy = m_K*(dist-l0)*(dist-l0);
    double m_Force = -m_K*(dist-l0);
}
void SpringPotentialBetweenTwoGroups::MovingVertex(vertex* v, Vec3D Dx)
{
    if(v->GetGroupName() == m_Group1Name)
    {
        m_Group1COG = m_Group1COG + Dx*(1/double(m_pGroup1.size()));
    }
    if(v->GetGroupName() == m_Group2Name)
    {
        m_Group2COG = m_Group2COG + Dx*(1/double(m_pGroup2.size()));
    }
}
void SpringPotentialBetweenTwoGroups::RejectMovingVertex(vertex* v, Vec3D Dx)
{
    if(v->GetGroupName() == m_Group1Name)
    {
        m_Group1COG = m_Group1COG - Dx*(1/double(m_pGroup1.size()));
    }
    if(v->GetGroupName() == m_Group2Name)
    {
        m_Group2COG = m_Group2COG - Dx*(1/double(m_pGroup2.size()));
    }
}
void SpringPotentialBetweenTwoGroups::MakeGroups(std::vector<vertex *> Ver,std::string ndx)
{

    
    m_Group1NDX = ReadIndex(ndx, m_Group1Name);
    m_Group2NDX = ReadIndex(ndx, m_Group2Name);
    
    for (int i=0;i<m_Group1NDX.size();i++)
    {
        if(m_Group1NDX[i]>Ver.size()-1)
        {
            std::cout<<"---->Error: number in the index does not match with the id of any vertex "<<std::endl;
            exit(0);
        }
        m_pGroup1.push_back(Ver.at(m_Group1NDX.at(i)));
        Ver.at(m_Group1NDX.at(i))->UpdateGroupName(m_Group1Name);
    }
    for (int i=0;i<m_Group2NDX.size();i++)
    {
        if(m_Group2NDX[i]>Ver.size()-1)
        {
            std::cout<<"---->Error: number in the index does not match with the id of any vertex "<<std::endl;
            exit(0);
        }
        m_pGroup2.push_back(Ver.at(m_Group2NDX.at(i)));
        Ver.at(m_Group2NDX.at(i))->UpdateGroupName(m_Group2Name);

    }
    m_Group1COG = COMVertexGroup(m_pGroup1);
    m_Group2COG = COMVertexGroup(m_pGroup2);
    Vec3D *pBox = (m_pGroup1.at(0))->GetBox();
    Vec3D PX1 = m_Group1COG;
    Vec3D PX2 = m_Group2COG;
    Vec3D Dist;
    
    for (int i=0;i<3;i++)
    {
        double dist=PX1(i)-PX2(i);
        if(dist>(*pBox)(i)/2.0)
        {
            dist=(*pBox)(i)-dist;
        }
        else if(dist<0)
        {
            dist=-dist;
        }
        Dist(i)=dist;
    }
    m_L0=m_Direction(0)*Dist(0)*Dist(0)+m_Direction(1)*Dist(1)*Dist(1)+m_Direction(2)*Dist(2)*Dist(2);
    m_L0=sqrt(m_L0);
}
/// Funcation to read index file
std::vector<int> SpringPotentialBetweenTwoGroups::ReadIndex(std::string ndx, std::string groupname)
{
    Nfunction f;
    std::string ext = ndx.substr(ndx.find_last_of(".") + 1);
    std::string filename;
    
    if(ndx.at(ndx.size()-1)=='x' && ndx.at(ndx.size()-2)=='n' && ndx.at(ndx.size()-3)=='i' && ndx.at(ndx.size()-4)=='.' )
        filename = ndx;
    else
        filename = ndx+".inx";
    
    if(f.FileExist(filename)==false)
    {
        std::cout<<"-----> Error: the index file name with the name "<<filename<<" does not exist"<<std::endl;
        exit(0);
    }
    std::ifstream indexfile;
    indexfile.open(filename.c_str());

    std::vector<int> Gndx;
    bool groupexist = false;
    int NAtom = 0;
    while (true)
    {
        std::string w1;
        indexfile>>w1;
         if(w1==groupname)
          {
              groupexist = true;
              indexfile>>NAtom;
              for (int i=0;i<NAtom;i++)
              {
                  int id;
                  indexfile>>id;
                  Gndx.push_back(id);
              }
              break;
          }
         else
         {
            std::getline(indexfile,w1);
         }
        
        if(indexfile.eof())
            break;

    }
    if(groupexist==false)
    {
        std::cout<<"----> Error: group with the name "<<groupname<<" does not exist in the index file "<<std::endl;
        exit(0);
    }
    indexfile.close();
    return Gndx;
}
Vec3D SpringPotentialBetweenTwoGroups::COMVertexGroup(std::vector<vertex *> ver)
{
    Vec3D com;
    double x=0;
    double y=0;
    double z=0;
    for (std::vector<vertex *>::iterator it = ver.begin() ; it != ver.end(); ++it)
    {
        x+=(*it)->GetVXPos();
        y+=(*it)->GetVYPos();
        z+=(*it)->GetVZPos();
    }
    com(0)=x/ver.size();
    com(1)=y/ver.size();
    com(2)=z/ver.size();

    return com;
}

