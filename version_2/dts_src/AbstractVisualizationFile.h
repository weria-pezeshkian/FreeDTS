#if !defined(AFX_AbstractVisualizationFile_H)
#define AFX_AbstractVisualizationFile_H

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
========================================================
*/
class  AbstractVisualizationFile {
public:
    AbstractVisualizationFile(){
        m_Period = 10;
    }
    virtual ~ AbstractVisualizationFile(){
        
    }
    virtual bool WriteAFrame(int step) = 0;
    virtual bool OpenFolder() = 0;
    virtual inline  std::string GetDerivedDefaultReadName()  {return "";}
    virtual std::string CurrentState() = 0;
    inline int GetPeriod()  {return m_Period;}
    inline static std::string GetBaseDefaultReadName()  {return "VisualizationFormat";}
    

    void SetPeriod(int period){
        m_Period = period;
        return;
    }
    
protected:
    int m_Period;
    
private:
    
};
#endif
