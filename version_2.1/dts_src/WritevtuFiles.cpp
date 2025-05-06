#include <iostream>
#include <string.h>
#include "State.h"
#include "WritevtuFiles.h"

WritevtuFiles::WritevtuFiles(State* pState, int period, std::string foldername) :
    m_FolderName(foldername),
    m_Box(pState->GetMesh()->Link2ReferenceBox())
 {
    m_pState = pState;
    m_Period = period;
}
WritevtuFiles::WritevtuFiles(State* pState) :
    m_Box(pState->GetMesh()->Link2ReferenceBox())
{
    
    m_pState = pState;
    m_FolderName = "VTU_Frames";
    m_Period = 10;
}
WritevtuFiles::~WritevtuFiles(){
    
}
bool WritevtuFiles::OpenFolder(){
  
  //  m_pBox = (pState->GetMesh())->GetBox();
    return Nfunction::OpenFolder(m_FolderName);
}
void WritevtuFiles::WriteInclusion(std::string id, const std::vector<vertex *>  &all_ver, std::ofstream *Output)
{
    double cc=0;
    (*Output)<<std::fixed;
    (*Output)<<std::setprecision( Precision );
    *Output<<"        <DataArray type=\"Float32\" Name=\""<<id<<"\"  NumberOfComponents=\"3\" Format=\"ascii\">"<<"\n";
    
   // std::vector<vertex *>  all_ver = m_pState->GetMesh()->GetActiveV();
    for (std::vector<vertex*>::const_iterator itr = all_ver.begin() ; itr != all_ver.end(); ++itr)
    {
        if((*itr)->VertexOwnInclusion()){
            
            Vec3D direction = (*itr)->GetInclusion()->GetLDirection();
            direction=((*itr)->GetL2GTransferMatrix())*direction;
            
            *Output<<"  "<<direction(0)<<"      "<<direction(1)<<"      "<<direction(2)<<"\n";
        }
        else
        *Output<<"  "<<cc<<"    "<<cc<<"    "<<cc<<"\n";
    }
    *Output<<"        </DataArray>"<<"\n";
    
}

bool WritevtuFiles::WriteVectorFields(const std::vector<vertex *> &all_ver, std::ofstream *Output) {
    /*
     * @brief Write vector fields to a VTU file.
     *
     * This function writes the vector fields of all vertices to a VTU file.
     *
     * @param all_ver A vector containing pointers to all vertex objects.
     * @param Output A pointer to the output file stream.
     * @return True if the operation is successful, false otherwise.
     */
    // Set output format to fixed and precision to specified value
    
    int number_of_vector_fields = m_pState->GetMesh()->GetNoVFPerVertex();
    if(number_of_vector_fields == 0){
        return true;
    }
    
    (*Output) << std::fixed;
    (*Output) << std::setprecision(Precision);

    // Loop through each vector field per vertex
    for (int i = 0; i < number_of_vector_fields; i++) {
        // Create the name ID for the vector field
        std::string nameid = "Vector_Field" + Nfunction::D2S(i + 1);

        // Write the DataArray header to the output file
        *Output << "        <DataArray type=\"Float32\" Name=\"" << nameid
                << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

        // Loop through all vertices and write their vector field data
        for (std::vector<vertex*>::const_iterator itr = all_ver.begin(); itr != all_ver.end(); ++itr) {
            // Get the local direction of the vector field and transform it to global coordinates
            Vec3D direction = (*itr)->GetVectorField(i)->GetLDirection();
            direction = ((*itr)->GetL2GTransferMatrix()) * direction;

            // Write the transformed vector field data to the output file
            *Output << "  " << direction(0) << "      " << direction(1) << "      " << direction(2) << "\n";
        }

        // Close the DataArray element
        *Output << "        </DataArray>\n";
    }

    return true;
}
bool WritevtuFiles::WriteAFrame(int step){
    

    if(m_Period == 0 || step%m_Period != 0){
        return false;
    }
    int file_index = step/m_Period;
    std::string Filename = "./"+m_FolderName+"/conf"+Nfunction::Int_to_String(file_index)+".vtu";
    
    std::ofstream Output;
    Output.open(Filename.c_str());
    
    //m_pBox = m_pState->GetMesh()->GetBox();
    double Half_Lx = m_Box(0)/2.0;
    double Half_Ly = m_Box(1)/2.0;
    double Half_Lz = m_Box(2)/2.0;

    
    std::vector<triangle *>  all_tri = m_pState->GetMesh()->GetActiveT();
    std::vector<vertex *>  all_ver = m_pState->GetMesh()->GetActiveV();
    int numv=all_ver.size();
    int numtri=all_tri.size();
    int numtrirep=0;
    
    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it){
        
        (*it)->UpdateRepresentation(true);
        double dx1=(*it)->GetV1()->GetXPos()-(*it)->GetV2()->GetXPos();
        double dy1=(*it)->GetV1()->GetYPos()-(*it)->GetV2()->GetYPos();
        double dz1=(*it)->GetV1()->GetZPos()-(*it)->GetV2()->GetZPos();
        double dx2=(*it)->GetV1()->GetXPos()-(*it)->GetV3()->GetXPos();
        double dy2=(*it)->GetV1()->GetYPos()-(*it)->GetV3()->GetYPos();
        double dz2=(*it)->GetV1()->GetZPos()-(*it)->GetV3()->GetZPos();
        
        if(fabs(dx1)>Half_Lx || fabs(dy1)>Half_Ly || fabs(dz1)>Half_Lz){
            (*it)->UpdateRepresentation(false);
        }
        else if(fabs(dx2)>Half_Lx || fabs(dy2)>Half_Ly || fabs(dz2)>Half_Lz){
                (*it)->UpdateRepresentation(false);
        }
        else{
            numtrirep++;
        }
    }
    Output<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">"<<"\n";
    Output<<"  <UnstructuredGrid>"<<"\n";
    Output<<"    <Piece NumberOfPoints=\""<<numv<<"\" NumberOfCells=\""<<numtrirep<<"\">"<<"\n";
    Output<<"      <PointData Scalars=\"scalars\">"<<"\n";
        
        

    std::vector <InclusionType*> TYPES = (m_pState->GetMesh())->GetInclusionType();
    std::vector<int> inctypeid;
    std::vector<std::string> inctype;
    for (std::vector<InclusionType*>::iterator it = TYPES.begin() ; it != TYPES.end(); ++it){
        
        inctypeid.push_back((*it)->ITid);
        inctype.push_back((*it)->ITName);
    }

        
    for (int n=0;n<inctypeid.size();n++) {
    Output<<"        <DataArray type=\"Float32\" Name=\""<<inctype[n]<<"\" Format=\"ascii\">"<<"\n";
        for (std::vector<vertex *>::iterator it = all_ver.begin() ; it != all_ver.end(); ++it){
        
            if((*it)->VertexOwnInclusion()==true){
                if((*it)->GetInclusion()->GetInclusionType()->ITid ==inctypeid[n]){
                    Output<<"          "<<1<<"\n";
                }
                else{
                    Output<<"          "<<0<<"\n";
                }
            }
            else{
                Output<<"          "<<0<<"\n";
            }

        }
        Output<<"        </DataArray>"<<"\n";
    }
    //=== write group id
    Output<<"        <DataArray type=\"Float32\" Name=\""<<"GroupID"<<"\" Format=\"ascii\">"<<"\n";
    for (std::vector<vertex *>::iterator it = all_ver.begin() ; it != all_ver.end(); ++it){
        Output<<"          "<<(*it)->GetGroup()<<"\n";

    }
     Output<<"        </DataArray>"<<"\n";
     WriteInclusion("dir", all_ver, &Output);
    
    WriteVectorFields(all_ver, &Output);
    
    for (int n=0; n<m_AllVectors.size();n++){
        WriteVector(m_AllVectorNames[n],m_AllVectors[n], &Output);
    }


    Output<<"      </PointData>"<<"\n";
    Output<<"      <Points>"<<"\n";
    Output<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">"<<"\n";

    for (std::vector<vertex*>::iterator it = all_ver.begin() ; it != all_ver.end(); ++it) {
    Output<<"          "<<(*it)->GetXPos()<<" "<<(*it)->GetYPos()<<" "<<(*it)->GetZPos()<<" "<<"\n";
    }
    Output<<"        </DataArray>"<<"\n";
    Output<<"      </Points>"<<"\n";
    Output<<"      <Cells>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">"<<"\n";


    for (std::vector<triangle *>::iterator it = all_tri.begin() ; it != all_tri.end(); ++it){

        if((*it)->GetRepresentation()){
        Output<<"           "<<(*it)->GetV1()->GetVID()<<" "<<(*it)->GetV2()->GetVID()<<" "<<(*it)->GetV3()->GetVID()<<" "<<"\n";
        }
    }
    
    Output<<"        </DataArray>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">"<<"\n";
    int ofset=3;
    Output<<"          ";
    for (int i=0;i<numtrirep;i++){

        Output<<ofset+3*i<<" ";
    }
    Output<<"\n";

    Output<<"        </DataArray>"<<"\n";
    Output<<"        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">"<<"\n";
    Output<<"          ";
    
    for (int i=0;i<numtrirep;i++){
        Output<<5<<" ";
    }
    Output<<"\n";

    Output<<"        </DataArray>"<<"\n";
    Output<<"      </Cells>"<<"\n";
    Output<<"    </Piece>"<<"\n";
    Output<<"  </UnstructuredGrid>"<<"\n";
    Output<<"</VTKFile> "<<"\n";
    
    Output.flush();
    Output.close();

    return true;
}
bool WritevtuFiles::WriteVector(const std::string &name, const std::vector<Vec3D >  &vecs, std::ofstream *Output){

    
    (*Output) << std::fixed;
    (*Output) << std::setprecision(Precision);

        // Write the DataArray header to the output file
        *Output << "        <DataArray type=\"Float32\" Name=\"" << name
                << "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";

        // Loop through all vertices and write their vector field data
        for (std::vector<Vec3D>::const_iterator itr = vecs.begin(); itr != vecs.end(); ++itr) {
            *Output << "  " << (*itr)(0) << "      " << (*itr)(1) << "      " << (*itr)(2) << "\n";
        }

        // Close the DataArray element
        *Output << "        </DataArray>\n";


    return true;
    
}
bool WritevtuFiles::AddVector(const std::string& name, const std::vector<Vec3D >  &vecs){
    m_AllVectors.push_back(vecs);
    m_AllVectorNames.push_back(name);
    
    return true;
}
bool WritevtuFiles::ClearVector(){
    m_AllVectors.clear();
    m_AllVectorNames.clear();
    
    return true;
}
std::string WritevtuFiles::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName() +" "+m_FolderName+" "+ Nfunction::D2S(m_Period);
    return state;
}
