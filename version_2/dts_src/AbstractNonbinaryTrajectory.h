#if !defined(AFX_AbstractNonbinaryTrajectory_H)
#define AFX_AbstractNonbinaryTrajectory_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
========================================================
*/
class State;
class  AbstractNonbinaryTrajectory {
public:
    AbstractNonbinaryTrajectory(){
        
    }
    virtual ~ AbstractNonbinaryTrajectory(){
        
    }
    virtual void WriteAFrame(int step) = 0;
    virtual bool OpenFolder()= 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;
    inline static std::string GetBaseDefaultReadName() {return "NonbinaryTrajectory";}
    
private:
    
};
#endif
