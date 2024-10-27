
#include "State.h"
#include "RigidWallTypes.h"

TwoFlatParallelWall::TwoFlatParallelWall(State* pState, double thickness,  char direction){
    
    m_pState = pState;
    m_HalfThickness = thickness/2;
    m_MidPlane = 0;
    if(direction=='X'){
        m_Element = 0;
        m_Direction = "X";
    }
    else if(direction=='Y'){
        m_Element = 1;
        m_Direction = "Y";
    }
    if(direction=='Z'){
        m_Element = 2;
        m_Direction = "Z";
    }
    else {
        *(m_pState->GetTimeSeriesLog()) << "---> such a direction for wall is wrong \n";
    }
}
TwoFlatParallelWall::~TwoFlatParallelWall(){
    
}
void TwoFlatParallelWall::Initialize() {
    // Initialize mid-plane and maximum distance variables
    m_MidPlane = 0;
    double dist_max = 0;

    // Get all surface vertices
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    size_t numVertices = pAllVertices.size();

    // Compute the mid-plane position
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        m_MidPlane += (*it)->GetPos()(m_Element);
    }
    m_MidPlane /= static_cast<double>(numVertices);

    // Find the maximum distance from the mid-plane
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        double dist = (*it)->GetPos()(m_Element) - m_MidPlane;
        if (fabs(dist) > dist_max) {
            dist_max = fabs(dist);
        }
    }

    // Adjust the half thickness if needed
    if (dist_max > m_HalfThickness) {
        double newHalfThickness = dist_max + 0.1;
        std::cout << "---> Note: The provided thickness for the ParallelWall boundary is too small. Increasing it to " << newHalfThickness << "\n";
        *(m_pState->GetTimeSeriesLog()) << "---> Note: The provided thickness for the ParallelWall boundary is too small. Increasing it to " << newHalfThickness << "\n";
        m_HalfThickness = newHalfThickness;
    }
}
bool TwoFlatParallelWall::MoveHappensWithinTheBoundary(double dx, double dy, double dz, vertex* p_ver){
    
    Vec3D dr(dx,dy,dz);
    dr=dr+p_ver->GetPos();
    
    if(dr(m_Element)>m_MidPlane + m_HalfThickness || dr(m_Element)>m_MidPlane - m_HalfThickness)
        return false;

    return true;
}
std::string TwoFlatParallelWall::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(2*m_HalfThickness) +" "+m_Direction;

    return state;
}

//===== EllipsoidalShell implementations
EllipsoidalShell::EllipsoidalShell(State* pState, double thickness, double r, double a, double b, double c)  :
            m_pState(pState),
            m_R(r),
            m_HalfThickness(thickness/2),
            m_A(a),
            m_B(b),
            m_C(c)
{

    
    
    *(m_pState->GetTimeSeriesLog()) << "---> this class is incomplete \n";
    
}
EllipsoidalShell::~EllipsoidalShell(){
    
}
void EllipsoidalShell::Initialize() {

    double size = m_A * m_A + m_C * m_C + m_B * m_B;
    size = sqrt(size);
    m_A = m_A*sqrt(3)/size;
    m_B = m_B*sqrt(3)/size;
    m_C = m_C*sqrt(3)/size;

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    double V_number = double(pAllVertices.size());
    if(V_number == 0){
        std::cout<<"--> error: should not happen 22 \n";
    }
    V_number = 1.0/V_number;
     for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
         m_COG = m_COG + (*it)->GetPos()*(V_number);
      }
    
    double rmax = 0; double rmin = std::numeric_limits<double>::max();
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        
        Vec3D P = (*it)->GetPos() - m_COG;
        P(0) = P(0)*(1/m_A);
        P(1) = P(1)*(1/m_B);
        P(2) = P(2)*(1/m_C);

        double r2 = Vec3D::dot(P,P);
        
        if( r2 > rmax ){
            rmax = r2;
        }
        if( r2 < rmin){
            rmin = r2;
        }
    }
    rmax = sqrt(rmax);
    rmin = sqrt(rmin);

    if(rmax > m_R + m_HalfThickness || rmin < m_R - m_HalfThickness){
        
        m_HalfThickness = rmax - rmin;
        m_HalfThickness = m_HalfThickness/2;
        m_R = rmax - m_HalfThickness + 0.2;
        
        std::cout << "---> Note: The provided thickness or R for the EllipsoidalShell boundary is too small " <<m_HalfThickness<<" "<<m_R<< "\n";
        *(m_pState->GetTimeSeriesLog()) << "---> Note: The provided thickness or R for the EllipsoidalShell boundary is too small " <<m_HalfThickness<<" "<<m_R<< "\n";
        

        
    }


}
bool EllipsoidalShell::MoveHappensWithinTheBoundary(double dx, double dy, double dz, vertex* p_ver){
    
    Vec3D P(dx,dy,dz);
    P = P + p_ver->GetPos() - m_COG;
    P(0) = P(0)*(1/m_A);
    P(1) = P(1)*(1/m_B);
    P(2) = P(2)*(1/m_C);
    double r = sqrt(Vec3D::dot(P,P));
    
    if(r > m_R + m_HalfThickness || r < m_R - m_HalfThickness){
        return false;
    }
    
    return true;
}
std::string EllipsoidalShell::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(2*m_HalfThickness) +" "+ Nfunction::D2S(m_R) +" "+ Nfunction::D2S(m_A) +" "+ Nfunction::D2S(m_B) +" "+ Nfunction::D2S(m_C);
    return state;
}

//===== EllipsoidalCore implementations
EllipsoidalCore::EllipsoidalCore(State* pState, double r, double a, double b, double c)  :
            m_pState(pState),
            m_R(r),
            m_A(a),
            m_B(b),
            m_C(c)
{

    
    
}
EllipsoidalCore::~EllipsoidalCore(){
    
}
void EllipsoidalCore::Initialize() {

    double size = m_A * m_A + m_C * m_C + m_B * m_B;
    size = sqrt(size);
    m_A = m_A*sqrt(3)/size;
    m_B = m_B*sqrt(3)/size;
    m_C = m_C*sqrt(3)/size;

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    double V_number = double(pAllVertices.size());
    if(V_number == 0){
        std::cout<<"--> error: number of vertices is zero, should not happen \n";
    }
    V_number = 1.0/V_number;
     for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
         m_COG = m_COG + (*it)->GetPos()*(V_number);
      }
    m_A = 1/m_A;    m_B = 1/m_B;    m_C = 1/m_C;
    double rmin = std::numeric_limits<double>::max();
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        
        Vec3D P = (*it)->GetPos() - m_COG;
        P(0) = P(0)*m_A;
        P(1) = P(1)*m_B;
        P(2) = P(2)*m_C;

        double r2 = Vec3D::dot(P,P);
        
        if( r2 < rmin){
            rmin = r2;
        }
    }
    rmin = sqrt(rmin);
    if(m_R > rmin){
        m_R = rmin - 0.001;
        *(m_pState->GetTimeSeriesLog()) << " Radius of the boundry core was reduced to "<<m_R<<" so no vertex is within the core \n";

    }
    m_R = m_R*m_R;

}
bool EllipsoidalCore::MoveHappensWithinTheBoundary(double dx, double dy, double dz, vertex* p_ver){
    
    Vec3D P(dx,dy,dz);
    P = P + p_ver->GetPos() - m_COG;
    P(0) = P(0)*m_A;
    P(1) = P(1)*m_B;
    P(2) = P(2)*m_C;
    
    if(Vec3D::dot(P,P) < m_R ){
        return false;
    }
    
    return true;
}
std::string EllipsoidalCore::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(sqrt(m_R)) +" "+ Nfunction::D2S(1/m_A) +" "+ Nfunction::D2S(1/m_B) +" "+ Nfunction::D2S(1/m_C);
    return state;
}

