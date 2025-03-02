#if !defined(AFX_ActiveTwoStateInclusion_H_999B21B8_C13C_5648_BC23_444775086239__INCLUDED_)
#define AFX_ActiveTwoStateInclusion_H_999B21B8_C13C_5648_BC23_444775086239__INCLUDED_


#include "SimDef.h"
#include "AbstractInclusionConversion.h"
#include "inclusion.h"

/*
 * @class ActiveTwoStateInclusion
 * @brief Class for actively exchanging inclusion types between two states.
 *
 * This class implements an active inclusion exchange algorithm that continuously changes
 * the inclusion type between two specified states.
 *
 * @author Weria Pezeshkian
 */

class ActiveTwoStateInclusion : public AbstractInclusionConversion {

public:
    
    /*
     * @brief Constructor for ActiveTwoStateInclusion.
     * @param period The period of exchange cycles (on simulation steps).
     * @param ep .
     * @param per percentage of state one vs state two.
     * @param gama The activity parameter.
     * @param t_name1 The name of the first inclusion type.
     * @param t_name2 The name of the second inclusion type.
     */
    ActiveTwoStateInclusion(int period, double ep, double per, double gama, std::string t_name1,std::string t_name2 );
	 ~ActiveTwoStateInclusion();


    void Initialize(State *pstate);
    bool Exchange(int step);
    inline  std::string GetDerivedDefaultReadName()  {return "ActiveTwoStateInclusion";}
    inline static std::string GetDefaultReadName() {return "ActiveTwoStateInclusion";}
    std::string CurrentState();
    
    
private:
    int m_Period; ///< Period of exchange cycles.
    double m_ActiveEnergy; ///< total Active energy.
    double m_Epsilon; ///< epsilon is the rate, we could have different rate but reducing the number of the model parameters
    double m_Percentage; ///<  percentage of type 1 vs type 2 .
    double m_Gama; ///< Activity gamma parameter.
    int m_N2; ///< Number of inclusions of type 1.
    int m_N1; ///< Number of inclusions of type 2.
    int m_N; ///< Total number of inclusions.
    int m_Delta_N0; ///< Forced average of (n1 - n2).
    std::string m_TypeName_1; ///< Name of the first inclusion type.
    std::string m_TypeName_2; ///< Name of the second inclusion type.
    std::vector<inclusion*> m_pSubInc; ///< List of inclusions to be exchanged.
    InclusionType *m_pIncType1; ///< Pointer to the first inclusion type.
    InclusionType *m_pIncType2; ///< Pointer to the second inclusion type.
    ///
    ///
private:
    State* m_pState; ///< Pointer to the  state class.


};


#endif
