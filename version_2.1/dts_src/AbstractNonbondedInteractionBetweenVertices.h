#ifndef ABSTRACT_NONBONDED_INTERACTION_BETWEEN_VERTICES_H
#define ABSTRACT_NONBONDED_INTERACTION_BETWEEN_VERTICES_H

#include <iostream>
#include "vertex.h"

/*
=======================================================
 Developed 2024 by Weria Pezeshkian
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This is an abstract class for Nonbonded Interaction Between Vertices.

 Purpose:
 - Defines an abstract base class for nonbonded interactions.
 - Provides an interface for derived classes to implement specific interactions.
 - Includes a derived class that represents the absence of nonbonded interactions.

========================================================
*/

// Forward declaration of the State class
class State;

/**
 * @brief Abstract base class for nonbonded interactions between vertices.
 *
 * This class defines an interface for nonbonded interactions, requiring derived classes
 * to implement specific interaction behavior.
 */
class AbstractNonbondedInteractionBetweenVertices {
public:
    /**
     * @brief Default constructor.
     */
    AbstractNonbondedInteractionBetweenVertices() = default;

    /**
     * @brief Virtual destructor (ensures proper cleanup in derived classes).
     */
    virtual ~AbstractNonbondedInteractionBetweenVertices() = default;

    /**
     * @brief Computes the coupling energy for a given vertex.
     *
     * @param pvertex Pointer to a vertex object.
     * @return Computed coupling energy as a double.
     */
    virtual double GetVertexNonBondedEnergy(vertex* pvertex) const = 0;
    virtual void Initialize()  = 0;

    /**
     * @brief Retrieves the current state of the interaction.
     *
     * @return A string representation of the current state.
     */
    virtual std::string CurrentState() const = 0;

    /**
     * @brief Gets the derived class's default name for reading.
     *
     * @return A string representing the derived class's read name.
     */
    virtual std::string GetDerivedDefaultReadName() const = 0;

    /**
     * @brief Retrieves the base class's default read name.
     *
     * This function provides a default identifier for the base class.
     *
     * @return A string representing the base class's read name.
     */
    static inline std::string GetBaseDefaultReadName() { return "NonbondedInteractionBetweenVertices"; }
   // static inline std::string GetDefaultReadName() { return "NonbondedInteractionBetweenVertices"; }
    inline static std::string GetErrorMessage(std::string info) {
        return "---> error: unknown NonbondedInteractionBetweenVertices type -- \n";
    }

};

/**
 * @brief Derived class representing "No Nonbonded Interaction."
 *
 * This class implements the abstract interface but provides no actual nonbonded interaction.
 */
class NoNonbondedInteractionBetweenVertices : public AbstractNonbondedInteractionBetweenVertices {
public:
    /**
     * @brief Default constructor.
     */
    NoNonbondedInteractionBetweenVertices() = default;

    /**
     * @brief Destructor.
     */
    ~NoNonbondedInteractionBetweenVertices() override = default;

    inline static std::string GetDefaultReadName() {return "No";}
    void Initialize() override {return;};

    /**
     * @brief Computes the coupling energy for a given vertex.
     *
     * Since this class represents the absence of interaction, it always returns 0.
     *
     * @param pvertex Pointer to a vertex object.
     * @return Always returns 0.0.
     */
    double GetVertexNonBondedEnergy(vertex* pvertex) const override { return 0.0; }

    /**
     * @brief Retrieves the current state of the interaction.
     *
     * Constructs a string representation of the interaction state.
     *
     * @return A string in the format "NonbondedInteractionBetweenVertices = No".
     */
    std::string CurrentState() const override {
        std::string state_here = this->GetBaseDefaultReadName() + " = " + GetDerivedDefaultReadName();
        return state_here;
    }

    /**
     * @brief Gets the default read name for this derived class.
     *
     * @return Always returns "No".
     */
    std::string GetDerivedDefaultReadName() const override { return "No"; }
};

#endif // ABSTRACT_NONBONDED_INTERACTION_BETWEEN_VERTICES_H
