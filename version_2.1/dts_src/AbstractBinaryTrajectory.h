#if !defined(AFX_AbstractBinaryTrajectory_H)
#define AFX_AbstractBinaryTrajectory_H
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
class  AbstractBinaryTrajectory {
public:
    AbstractBinaryTrajectory(){
        
    }
    virtual ~ AbstractBinaryTrajectory(){
        
    }
    virtual bool OpenFile(bool clear, char rw) = 0;
    virtual void WriteAFrame(int &step) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "BinaryTrajectory";}

    
private:
    
};
class NoFile : public AbstractBinaryTrajectory {
public:
    NoFile(){
        
    }
    ~NoFile(){
        
    }
    
    virtual inline std::string GetDerivedDefaultReadName() {return "NoFile";}
    inline bool OpenFile(bool clear, char rw){
        return true;
    }

    
    void WriteAFrame(int &step){
        return;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};
#endif
