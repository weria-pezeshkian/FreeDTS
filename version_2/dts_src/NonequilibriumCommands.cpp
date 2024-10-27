#include "NonequilibriumCommands.h"
#include "State.h"
#include "RigidWallTypes.h"
#include <iostream>
#include <sstream>

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

    return true;
}

void NonequilibriumCommands::ExpandEllipsoidalCoreWall(int rate, double dr) {
    
    if(m_ActiveSimStep%rate != 0){
        return;
    }
    AbstractBoundary* pBoundary = m_pState->GetBoundary();
    if (EllipsoidalCore* pB = dynamic_cast<EllipsoidalCore*>(pBoundary)) {
        pB->m_R = pB->m_R + dr;
    }
    else {
        std::cout << "---> warning (maybe an error): you have called changes in EllipsoidalCoreWall while such a wall has not been called." << std::endl;
    }
}
