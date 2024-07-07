

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Convert.h"
#include "TSIFile.h"
#include "QFile.h"
#include "WritevtuFiles.h"
#include "VMDOutput.h"

Convert::Convert()
{
}
Convert::Convert(std::vector <std::string> argument)
{

    m_Argument = argument;
    m_Healthy =true;
    m_Seed =36723;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_OutputFilename = "out.q";
    m_InputFilename = "in.tsi";
    m_Box(0) = 10;
    m_Box(1) = 10;
    m_Box(2) = 10;
    m_Zoom(0) = 1;
    m_Zoom(1) = 1;
    m_Zoom(2) = 1;
    m_center = false;
    m_Translate(0)=0; m_Translate(1)=0; m_Translate(2)=0;
    Nfunction f;
    ExploreArguments();     // read the input data

    std::string ext = m_InputFilename.substr(m_InputFilename.find_last_of(".") + 1);
    QFile Q;
    TSIFile tsi;
    WritevtuFiles vtu;
    VMDOutput  gro;
    if(ext==TSExt)
    {
        m_BluePrint = Q.Read(m_InputFilename);
    }
    else if(ext==TSIExt)
    {
        m_BluePrint = tsi.ReadTSI(m_InputFilename);
    }
    else
    {
        std::cout<<"---> Error: input file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }
    //==============================================================
    
    // make the changes into the mesh
    
        std::vector<Vertex_Map> bvertex = m_BluePrint.bvertex;
        m_Box = m_BluePrint.simbox;
        ExploreArguments();     // read the input data
    
    // first we do the translation
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
         double x = it->x + m_Translate(0);
         double y = it->y + m_Translate(1);
         double z = it->z + m_Translate(2);
            x=x*(m_Zoom(0));
            y=y*(m_Zoom(1));
            z=z*(m_Zoom(2));
        m_Box(0) = m_Box(0)*fabs(m_Zoom(0));
        m_Box(1) = m_Box(1)*fabs(m_Zoom(1));
        m_Box(2) = m_Box(2)*fabs(m_Zoom(2));

         if(x>m_Box(0))
             x = x-(int(x/m_Box(0)))*m_Box(0);
         else if(x<0)
             x= x+(int(fabs(x)/m_Box(0))+1)*m_Box(0);
        
        if(y>m_Box(1))
            y = y-(int(y/m_Box(1)))*m_Box(1);
        else if(y<0)
            y= y+(int(fabs(y)/m_Box(1))+1)*m_Box(1);
        
        if(z>m_Box(2))
            z = z-(int(z/m_Box(2)))*m_Box(2);
        else if(z<0)
            z= z+(int(fabs(z)/m_Box(2))+1)*m_Box(2);
        
        it->x = x;
        it->y = y;
        it->z = z;
        

    }
    // we do box change
    m_BluePrint.simbox = m_Box;
    
    // we do centering;
    double xcm = 0;    double ycm = 0; double zcm = 0;

    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        xcm+= it->x;
        ycm+= it->y;
        zcm+= it->z;
    }
    xcm = xcm/double(bvertex.size());
    ycm = ycm/double(bvertex.size());
    zcm = zcm/double(bvertex.size());

    if(m_center == true)
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        it->x = it->x - xcm + m_Box(0)/2;
        it->y = it->y - ycm + m_Box(1)/2;
        it->z = it->z - zcm + m_Box(2)/2;
        
    }
    m_BluePrint.bvertex = bvertex;

    

    
    //=======================================================
    std::string outext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    if(outext==TSExt)
    {
        Q.Write(m_OutputFilename, m_BluePrint);
    }
    else if(outext==TSIExt)
    {
        tsi.WriteTSI(m_OutputFilename, m_BluePrint);
    }
    else if(outext=="vtu")
    {
        vtu.Write(m_BluePrint, m_OutputFilename);
    }
    else if(outext=="gro")
    {
        gro.Write(m_OutputFilename, m_BluePrint);
    }
    else
    {
        std::cout<<"---> Error: input file with "<<ext<<" extension is not recognized. "<<std::endl;
    }

    
}
Convert::~Convert()
{
    
}
void Convert::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-o")
        {
            m_OutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-in")
        {
            m_InputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-zoom")
        {
            m_Zoom(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Zoom(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Zoom(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-translate")
        {
            m_Translate(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Translate(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Translate(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;

        }
        else if(Arg1=="-Box")
        {
            m_Box(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Box(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Box(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
            
        }
        else if(Arg1=="-center")
        {
            m_center = true;
            i=i-1;
            
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else
        {
            std::cout << "---> Error: "<<Arg1<<" ... bad argument ";
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }

}
void Convert::HelpMessage()
{
    std::cout<<"--------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"--------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting  -------------------"<<"\n";
    std::cout<<" ./CNV -in out.q  -o out.tsi"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -in         string       out.q               input file name, could be tsi or q file formats "<<"\n";
    std::cout<<"  -o          string       out.tsi             output file name, could be tsi/vtu/gro or q file formats "<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
}
