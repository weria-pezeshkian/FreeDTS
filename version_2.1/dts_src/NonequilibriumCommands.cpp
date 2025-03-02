#include <iostream>
#include <sstream>
#include "NonequilibriumCommands.h"
#include "State.h"
#include "RigidWallTypes.h"


NonequilibriumCommands::NonequilibriumCommands(State* pState) {
    m_pState = pState;
    m_ActiveSimStep = 0;
}

NonequilibriumCommands::~NonequilibriumCommands() {
}

void NonequilibriumCommands::Run(int step) {
    m_ActiveSimStep = step;

    // Iterate over the container and call each lambda
    for (std::size_t i = 0; i < m_FunctionContainer.size(); ++i) {
        // Execute each stored function (with captured arguments)
        m_FunctionContainer[i]();
    }
}

bool NonequilibriumCommands::LoadCommand(std::string strcommand) {
    std::vector<std::string> command_arguments = Nfunction::Split(strcommand);
    if (command_arguments.size() < 2) {
        std::cout << "  error-> NonequilibriumCommands: not enough arguments in the input file" << std::endl;
        return false;
    }

    if (command_arguments[0] == "ExpandEllipsoidalCoreWall") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 3) {
            std::cout << "  error-> NonequilibriumCommands ExpandEllipsoidalCoreWall: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double Dr = Nfunction::String_to_Double(command_arguments[2]);

        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, Dr]() { ExpandEllipsoidalCoreWall(Rate, Dr); });
    }
    else if (command_arguments[0] == "ThinningEllipsoidalShell") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 3) {
            std::cout << "  error-> NonequilibriumCommands ExpandEllipsoidalCoreWall: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double Dr = Nfunction::String_to_Double(command_arguments[2]);

        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, Dr]() { ThinningEllipsoidalShell(Rate, Dr); });
    }
    else if (command_arguments[0] == "IncrementHarmonicPotentialBetweenTwoGroups") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 3) {
            std::cout << "  error-> NonequilibriumCommands IncrementHarmonicPotentialBetweenTwoGroups: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double Dr = Nfunction::String_to_Double(command_arguments[2]);

        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, Dr]() { IncrementHarmonicPotentialBetweenTwoGroups(Rate, Dr); });
    }
    else if (command_arguments[0] == "IncrementSphericalSubstrateCenter") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 6) {
            std::cout << "  error-> NonequilibriumCommands IncrementHarmonicPotentialBetweenTwoGroups: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double Dr = Nfunction::String_to_Double(command_arguments[2]);
        double nx = Nfunction::String_to_Double(command_arguments[3]);
        double ny = Nfunction::String_to_Double(command_arguments[4]);
        double nz = Nfunction::String_to_Double(command_arguments[5]);
        Vec3D Dirct(nx,ny,nz);
        Dirct.normalize();
        Dirct*Dr;
        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, Dirct]() { IncrementSphericalSubstrateCenter(Rate, Dirct); });
    }
    else if (command_arguments[0] == "IncrementVolumeCouplingSecondOrder") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 3) {
            std::cout << "  error-> NonequilibriumCommands IncrementVolumeCouplingSecondOrder: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double Dr = Nfunction::String_to_Double(command_arguments[2]);

        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, Dr]() { IncrementVolumeCouplingSecondOrder(Rate, Dr); });
    }
    else if (command_arguments[0] == "ChangeTemperatureConstantRate") {
        // Convert the second argument to an integer (e.g., X value)
        if (command_arguments.size() < 5) {
            std::cout << "  error-> NonequilibriumCommands ChangeTemperatureConstantRate: not enough arguments in the input file" << std::endl;
            return false;
        }

        int Rate = Nfunction::String_to_Int(command_arguments[1]);
        double DT = Nfunction::String_to_Double(command_arguments[2]);

        // Use lambda to store the function and argument
        m_FunctionContainer.push_back([this, Rate, DT]() { ChangeTemperatureWithConstantRate(Rate, DT ); });
    }

    return true;
}
void NonequilibriumCommands::ThinningEllipsoidalShell(int rate, double dr) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }

    // Get boundary from state
    AbstractBoundary* pBoundary = m_pState->GetBoundary();
    
    // Check if the boundary is an EllipsoidalCore
    if (EllipsoidalShell* pB = dynamic_cast<EllipsoidalShell*>(pBoundary)) {
        
        // Check if all vertices are within the boundary
        bool allVerticesInside = std::all_of(m_pState->GetMesh()->GetActiveV().begin(),
                                             m_pState->GetMesh()->GetActiveV().end(),
                                             [&pB](vertex* v) {
                                                 return pB->MoveHappensWithinTheBoundary(0, 0, 0, v);
                                             });

        // If all vertices are inside the boundary, update the radius
        if (allVerticesInside) {
            pB->m_HalfThickness -= dr;

        }

    } else {
        // Error message if the boundary is not an EllipsoidalCore
        std::cerr << "---> Error: Attempted to expand an EllipsoidalShell when no EllipsoidalShell is present in the boundary." << std::endl;
    }
}
void NonequilibriumCommands::ExpandEllipsoidalCoreWall(int rate, double dr) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }

    // Get boundary from state
    AbstractBoundary* pBoundary = m_pState->GetBoundary();
    
    // Check if the boundary is an EllipsoidalCore
    if (EllipsoidalCore* pB = dynamic_cast<EllipsoidalCore*>(pBoundary)) {
        
        // Check if all vertices are within the boundary
        bool allVerticesInside = std::all_of(m_pState->GetMesh()->GetActiveV().begin(),
                                             m_pState->GetMesh()->GetActiveV().end(),
                                             [&pB](vertex* v) {
                                                 return pB->MoveHappensWithinTheBoundary(0, 0, 0, v);
                                             });

        // If all vertices are inside the boundary, update the radius
        if (allVerticesInside) {
            pB->m_R += dr;
        }

    } else {
        // Error message if the boundary is not an EllipsoidalCore
        std::cerr << "---> Error: Attempted to expand an EllipsoidalCoreWall when no EllipsoidalCore is present in the boundary." << std::endl;
    }
}
void NonequilibriumCommands::IncrementHarmonicPotentialBetweenTwoGroups(int rate, double dr) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }

    // Get boundary from state
    AbstractApplyConstraintBetweenGroups* pConstraintBetweenGroups = m_pState->GetApplyConstraintBetweenGroups();
    
    // Check if the boundary is an EllipsoidalCore
    if (HarmonicPotentialBetweenTwoGroups* pB = dynamic_cast<HarmonicPotentialBetweenTwoGroups*>(pConstraintBetweenGroups)) {
            pB->m_L0 += dr;

    } else {
        // Error message if the boundary is not an EllipsoidalCore
        std::cerr << "---> Error: Attempted to expand an IncrementHarmonicPotentialBetweenTwoGroups when not applied." << std::endl;
    }
}
void NonequilibriumCommands::IncrementVolumeCouplingSecondOrder(int rate, double dr) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }

    // Get boundary from state
    AbstractVolumeCoupling* pVolumeCouplingSecondOrder = m_pState->GetVolumeCoupling();
    
    // Check if the boundary is an EllipsoidalCore
    if (VolumeCouplingSecondOrder* pB = dynamic_cast<VolumeCouplingSecondOrder*>(pVolumeCouplingSecondOrder)) {
            pB->m_TargetV += dr;

    } else {
        // Error message if the boundary is not an EllipsoidalCore
        std::cerr << "---> Error: Attempted to expand an IncrementHarmonicPotentialBetweenTwoGroups when not applied." << std::endl;
    }
}
void NonequilibriumCommands::IncrementSphericalSubstrateCenter(int rate, Vec3D Dirct) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }

    // Get boundary from state
    AbstractVertexAdhesionToSubstrate* pVertexAdhesionToSubstrate = m_pState->GetVertexAdhesionToSubstrate();
    
    // Check if the boundary is an EllipsoidalCore
    if (SphericalVertexSubstrate* pB = dynamic_cast<SphericalVertexSubstrate*>(pVertexAdhesionToSubstrate)) {
        pB->m_Center = pB->m_Center + Dirct;

    } else {
        // Error message if the boundary is not an EllipsoidalCore
        std::cerr << "---> Error: Attempted to expand an IncrementHarmonicPotentialBetweenTwoGroups when not applied." << std::endl;
    }
}
void NonequilibriumCommands::ChangeTemperatureWithConstantRate(int rate, double DT) {

    // Skip steps based on rate
    if (m_ActiveSimStep % rate != 0) {
        return;
    }
    double beta =  m_pState->GetSimulation()->GetBeta();
    beta += DT;
    m_pState->GetSimulation()->SetBeta(beta, m_pState->GetSimulation()->GetDBeta());

    //    double m_MinLength2;
    //double m_MaxLength2;
}
