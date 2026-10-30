#ifndef EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H
#define EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H

#include <string>
#include <vector>

#include "SimDef.h"
#include "AbstractInclusionConversion.h"
#include "inclusion.h"
//class State;

/**
 * @class EquilibriumInclusionExchangeByChemicalPotential
 * @brief Performs inclusion exchanges between two types based on chemical potential equilibrium.
 *
 * This class implements an algorithm that exchanges inclusions between two
 * specified types within a simulation mesh according to a target chemical potential.
 * The exchanges are probabilistic and guided by the energetics of the system,
 * ensuring that the process respects equilibrium thermodynamics.
 *
 * The class handles:
 *  - Initialization of inclusion lists of the target types.
 *  - Calculation of the number of attempted moves per exchange cycle.
 *  - Evaluation of the energy change for each potential inclusion swap.
 *  - Probabilistic acceptance or rejection of each swap using a Monte Carlo scheme.
 *  - Tracking the current number of inclusions of each type.
 *
 * Exchanges are performed at user-defined intervals (`m_Period`) and with
 * user-defined intensity (`m_Rate`). Chemical potential differences (`m_Mu`)
 * drive the exchange process to reach equilibrium.
 *
 * @note The algorithm assumes that only two inclusion types are involved.
 *       Attempting to initialize with more than two types will result in undefined behavior.
 *
 * @author
 * Weria Pezeshkian
 */

class EquilibriumInclusionExchangeByChemicalPotential : public AbstractInclusionConversion
{
public:
    /**
     * @brief Constructor
     * @param pstate Pointer to the current simulation state
     * @param period Number of simulation steps between successive exchange attempts
     * @param rate Fraction of inclusions to attempt swapping per cycle
     * @param mu Chemical potential difference between the two inclusion types
     * @param type1 Name of the first inclusion type
     * @param type2 Name of the second inclusion type
     *
     * Initializes the object with simulation parameters and target inclusion types.
     * Sets up counters and links to the state-dependent thermodynamic parameters.
     */
    EquilibriumInclusionExchangeByChemicalPotential(int period,
                                                    double rate,
                                                    double mu,
                                                    std::string type1,
                                                    std::string type2);

    ~EquilibriumInclusionExchangeByChemicalPotential() override;

    void Initialize() override;
    bool Exchange(int step) override;

    std::string CurrentState() override;

    static inline std::string GetDefaultReadName() { return "EquilibriumExchange"; }
    inline std::string GetDerivedDefaultReadName() override {
        return GetDefaultReadName();
    }

    /**
     * @brief Attempts a single inclusion swap based on the Monte Carlo acceptance criterion.
     * @param p_inc Pointer to the inclusion to attempt swapping
     * @param thermal Random number between 0 and 1 for probabilistic acceptance
     * @return true if the move was accepted, false if rejected
     *
     * This method computes the energy change of swapping the inclusion type and
     * uses the chemical potential difference to determine the acceptance probability.
     * If accepted, the system energy is updated; otherwise, the inclusion is reverted.
     */
    bool TryForOneInclusion(inclusion * pinc, double thermal);

private:
    int m_Period;                   ///< Exchange cycle period
    double m_Mu;                    ///< Chemical potential
    double m_Rate;                    ///< Chemical potential

    int m_N1;                       ///< Number of type1 inclusions
    int m_N2;                       ///< Number of type2 inclusions
    int m_N;                        ///< Total number of inclusions

    int m_NumberOfMoves;
    std::string m_TypeName_1;        ///< Type name 1
    std::string m_TypeName_2;        ///< Type name 2

    std::vector<inclusion*> m_pSubInc; ///< Inclusions available for exchange

    InclusionType* m_pIncType1;     ///< Pointer to type 1 definition
    InclusionType* m_pIncType2;     ///< Pointer to type 2 definition
    double *m_Beta;
    double *m_DBeta;
};

#endif // EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H
