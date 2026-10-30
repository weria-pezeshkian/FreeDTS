/*
================================================================
FactoryVertexAdhesionToSubstrateMethod
================================================================

Singleton factory for creating instances of classes derived from

    AbstractVertexAdhesionToSubstrate

at runtime using a string identifier.

----------------------------------------------------------------
Creator Requirements
----------------------------------------------------------------

Each derived class must provide a static creator function:

    static AbstractVertexAdhesionToSubstrate* Create(std::istream& input);

The creator function is responsible for:

1. Reading parameters from the input stream
2. Constructing the derived class
3. Returning a raw pointer

Example:

    static AbstractVertexAdhesionToSubstrate* Create(std::istream& input)
    {
        double k;
        input >> k;
        return new MyAdhesion(k);
    }

----------------------------------------------------------------
Registration
----------------------------------------------------------------

Each derived class registers itself with the factory:

    FactoryVertexAdhesionToSubstrateMethod::Instance()
        .Register("MyAdhesion", MyAdhesion::Create);

----------------------------------------------------------------
Usage
----------------------------------------------------------------

Objects are created using:

    AbstractVertexAdhesionToSubstrate* obj =
        FactoryVertexAdhesionToSubstrateMethod::Instance()
            .Create(name, input);

If the identifier is unknown, nullptr is returned.

================================================================
*/

#pragma once

#include <unordered_map>
#include <functional>
#include <string>
#include <istream>
#include <stdexcept>
#include <iostream>

class AbstractVertexAdhesionToSubstrate;

class FactoryVertexAdhesionToSubstrateMethod
{
public:

    // Creator function returns a raw pointer
    typedef AbstractVertexAdhesionToSubstrate* (*Creator)(std::istream&);

    /*
    ------------------------------------------------------------
    Singleton Instance
    ------------------------------------------------------------
    */
    static FactoryVertexAdhesionToSubstrateMethod& Instance()
    {
        static FactoryVertexAdhesionToSubstrateMethod instance;
        return instance;
    }

    /*
    ------------------------------------------------------------
    Register
    ------------------------------------------------------------
    Registers a creator function with a string identifier.

    Throws std::runtime_error if the identifier already exists.
    */
    void Register(const std::string& name, Creator creator)
    {
        std::unordered_map<std::string, Creator>::iterator it =
            m_Creators.find(name);

        if (it != m_Creators.end())
        {
            throw std::runtime_error(
                "FactoryVertexAdhesionToSubstrateMethod: "
                "duplicate registration of '" + name + "'"
            );
        }

        m_Creators[name] = creator;
    }

    /*
    ------------------------------------------------------------
    Create
    ------------------------------------------------------------
    Creates a derived adhesion method.

    Parameters:
        name  : identifier of the derived class
        input : stream containing constructor parameters

    Returns:
        AbstractVertexAdhesionToSubstrate* (raw pointer)
        nullptr if identifier is unknown
    */
    AbstractVertexAdhesionToSubstrate* Create(
        const std::string& name,
        std::istream& input)
    {
        std::unordered_map<std::string, Creator>::iterator it =
            m_Creators.find(name);

        if (it == m_Creators.end())
        {
            return nullptr;
        }

        return it->second(input);
    }

private:

    FactoryVertexAdhesionToSubstrateMethod() {}
    ~FactoryVertexAdhesionToSubstrateMethod() {}

    FactoryVertexAdhesionToSubstrateMethod(
        const FactoryVertexAdhesionToSubstrateMethod&) = delete;

    FactoryVertexAdhesionToSubstrateMethod& operator=(
        const FactoryVertexAdhesionToSubstrateMethod&) = delete;

private:

    std::unordered_map<std::string, Creator> m_Creators;
};
