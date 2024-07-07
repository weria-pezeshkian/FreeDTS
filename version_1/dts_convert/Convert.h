#if !defined(AFX_Convert_H_INCLUDED_)
#define AFX_Convert_H_INCLUDED_
#include "Def.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "MeshBluePrint.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 */

class Convert
{
public:
    
	Convert(std::vector <std::string> argument);
    Convert();
    ~Convert();

private:

public:
    std::vector <std::string> m_Argument;
    int m_Seed;       // seed for random number generator
    double m_MinFaceAngle;          //  minimum angle between the face (smaller will results in error), this is the value of the cos
    double m_MinVerticesDistanceSquare; //  minimum distance allowed between two vertices  (smaller will results in error)
    double m_MaxLinkLengthSquare;       //  maximum distance allowed between two nighbouring vertices  (larger will results in error)
    std::string m_OutputFilename; //  output TS file
    Vec3D m_Box;
    bool m_Healthy;
    std::string m_GeneralOutputFilename; //  a general file flag for specific run
    std::string m_InputFilename;
    Vec3D Box;
private:
    void ExploreArguments();         // updates variables based on the command line arguments
    void HelpMessage();              // writes a help message
    
private:
    MeshBluePrint m_BluePrint;

private:
    bool m_center;
    Vec3D m_Zoom;
    Vec3D m_Translate;
    
    



};

#endif
