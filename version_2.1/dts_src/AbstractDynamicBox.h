#if !defined(AFX_DynamicBox_H)
#define AFX_DynamicBox_H
#include <iostream>
/**
 * @file AbstractDynamicBox.h
 * @brief Defines the AbstractDynamicBox base class and the NoBoxChange derived class for box size management.
 *
 * This file contains the declaration of the AbstractDynamicBox class, which provides a base interface for
 * managing changes in the simulation box size, and the NoBoxChange class, which implements a no-change policy for the box size.
 *
 * =======================================================
 * Developed 2024 by Weria Pezeshkian
 * Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 * =======================================================
 *
 * The AbstractDynamicBox class defines a standard interface for dynamically changing the simulation box size.
 * It includes methods for changing the box size, initializing the box parameters, and obtaining the current state of the box.
 * The class also provides functionality for updating the change rate (DR) and calculating the acceptance rate of box size changes.
 *
 * The NoBoxChange class inherits from AbstractDynamicBox and implements a no-change policy, meaning the box size remains constant
 * throughout the simulation. This class overrides the virtual functions from AbstractDynamicBox to provide specific behavior
 * indicating that no changes to the box size are performed.
 *
 * Usage:
 * - Inherit from AbstractDynamicBox to create custom box size management policies.
 * - Use NoBoxChange for scenarios where the box size should remain unchanged.
 *
 * Example:
 * @code
 * NoBoxChange noChange;
 * noChange.Initialize();
 * bool changed = noChange.ChangeBoxSize(step);
 * std::string state = noChange.CurrentState();
 * @endcode
 */

class AbstractDynamicBox {
public:
    /*
     * @brief Default constructor initializes the change rate (DR).
     */
    AbstractDynamicBox() : m_DR(0.01), m_NumberOfAttemptedMoves(0), m_AcceptedMoves(0) {}

    /*
     * @brief Virtual destructor.
     */
    virtual ~AbstractDynamicBox() {}

    /*
     * @brief Pure virtual function to change the box size at a given step.
     *
     * @param step The current simulation step.
     * @return true if the box size was changed, false otherwise.
     */
    virtual bool ChangeBoxSize(int step) = 0;

    /*
     * @brief Pure virtual function to initialize the box parameters.
     */
    virtual void Initialize() = 0;

    /*
     * @brief Pure virtual function to get the current state of the box.
     *
     * @return A string representing the current state.
     */
    virtual std::string CurrentState() = 0;

    /*
     * @brief Pure virtual function to get the derived class's default read name.
     *
     * @return A string representing the default read name of the derived class.
     */
    virtual std::string GetDerivedDefaultReadName() = 0;

    /*
     * @brief Gets the base class's default read name.
     *
     * @return A string representing the default read name of the base class.
     */
    inline static std::string GetBaseDefaultReadName() { return "Dynamic_Box"; }

    /*
     * @brief Gets the current change rate (DR).
     *
     * @return The current change rate.
     */
    inline double GetDR() { return m_DR; }

    /*
     * @brief Updates the change rate (DR).
     *
     * @param dr The new change rate.
     */
    void UpdateDR(double dr) { m_DR = dr; }

    /*
     * @brief Calculates the acceptance rate of box size changes.
     *
     * @param reset Whether to reset the acceptance and attempted moves counters.
     * @return The acceptance rate.
     */
    double GetAcceptanceRate(bool reset) {
        if (m_NumberOfAttemptedMoves == 0) {
            return 0.0;
        }
        double rate = static_cast<double>(m_AcceptedMoves) / m_NumberOfAttemptedMoves;
        if (reset) {
            m_NumberOfAttemptedMoves = 0;
            m_AcceptedMoves = 0;
        }
        return rate;
    }

protected:
    double m_NumberOfAttemptedMoves;  ///< Number of attempted moves.
    double m_AcceptedMoves;           ///< Number of accepted moves.
    double m_DR;                      ///< Change rate.
};
/**
 * @class NoBoxChange
 * @brief Derived class that implements a no-change policy for the simulation box size.
 *
 * This class inherits from AbstractDynamicBox and overrides the virtual functions to indicate that the box size
 * should remain constant throughout the simulation.
 */
class NoBoxChange : public AbstractDynamicBox {
public:
    NoBoxChange(){
        
    }
    ~NoBoxChange(){
        
    }

    inline  std::string GetDerivedDefaultReadName()  {return "No";}
    inline  static std::string GetDefaultReadName()  {return "No";}

    void Initialize(){
        return;
    }
    bool ChangeBoxSize(int step){
        return false;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};

#endif
