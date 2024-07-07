#include <chrono>
#include "Three_Edge_Scission.h"
#include "State.h"
#include "MESH.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */

Three_Edge_Scission::Three_Edge_Scission(int period, State *pState) :
                AbstractSimulation(*(pState->GetSimulation())),
                MESH(*(pState->GetMesh())),
                m_pState(pState),
                m_Period(period) {

}
Three_Edge_Scission::~Three_Edge_Scission(){
    
}
void Three_Edge_Scission::Initialize() {

    std::cout<<"---> note: three_Edge_Scission initalized"<<std::endl;
}
bool Three_Edge_Scission::MCMove(int step) {
    
    if(m_Period==0 ||  step%m_Period!=0)
        return false;
    
   // m_Beta = m_pState->GetSimulation()->GetBeta();
    
  //  double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    //AddToTotalEnergy
    
    
    MCScissionMove(step);
    MCFussionMove(step);
    
    
    return true;
}
bool Three_Edge_Scission::MCFussionMove(int step){
    return false;
}
bool Three_Edge_Scission::MCScissionMove(int step)
{

    

    if(m_pGhostT.size()<4 || m_pGhostL.size()<4){
        std::cout<<" --->note: the number of the links and trinagles in repository is not enough, restart the simulations \n";
        exit(0);
    }

   // auto start = std::chrono::steady_clock::now();

    // finding the pair list for cutting the trinagles
    std::vector<pair_pot_triangle> pair_list  = FindPotentialTriangles();

    for (std::vector<pair_pot_triangle>::iterator it = pair_list.begin() ; it != pair_list.end(); ++it){// loop over all the possible one
        double enew = 0;
        double eold = 0;
//---- effected vertex energy
        eold+= ((it->PT1).pv1)->GetEnergy();
        eold+= ((it->PT1).pv2)->GetEnergy();
        eold+= ((it->PT1).pv3)->GetEnergy();
        eold+= ((it->PT2).pv1)->GetEnergy();
        eold+= ((it->PT2).pv2)->GetEnergy();
        eold+= ((it->PT2).pv3)->GetEnergy();
        std::vector <links *> Clinks = it->ConnectingLinks;
        std::vector <triangle *> Ctriangles = it->ConnectingTriangles;
//---- inclusion interaction energy in the energy change
        //-- getting all the effected edges for interaction energy
        {
        std::vector <links *> nl11 = ((it->PT1).pv1)->GetVLinkList();
        std::vector <links *> nl12 = ((it->PT1).pv2)->GetVLinkList();
        std::vector <links *> nl13 = ((it->PT1).pv3)->GetVLinkList();
        std::vector <links *> nl21 = ((it->PT2).pv1)->GetVLinkList();
        std::vector <links *> nl22 = ((it->PT2).pv2)->GetVLinkList();
        std::vector <links *> nl23 = ((it->PT2).pv3)->GetVLinkList();
        std::vector <links *> AllEffectedEdge = nl11;
            AllEffectedEdge.insert(AllEffectedEdge.end(), nl12.begin(), nl12.end());
            AllEffectedEdge.insert(AllEffectedEdge.end(), nl13.begin(), nl13.end());
            AllEffectedEdge.insert(AllEffectedEdge.end(), nl21.begin(), nl21.end());
            AllEffectedEdge.insert(AllEffectedEdge.end(), nl22.begin(), nl22.end());
            AllEffectedEdge.insert(AllEffectedEdge.end(), nl23.begin(), nl23.end());
            //-- remove mirorlinks
            for (std::vector<links *>::iterator it = Clinks.begin() ; it != Clinks.end(); ++it){
                RemoveFromLinkList(*it,AllEffectedEdge);
            }
            RemoveFromLinkList((it->PT1).pl1,AllEffectedEdge);
            RemoveFromLinkList((it->PT1).pl2,AllEffectedEdge);
            RemoveFromLinkList((it->PT1).pl3,AllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl1,AllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl2,AllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl3,AllEffectedEdge);
        //---- inclusion interaction energy;
         for (std::vector<links *>::iterator it = AllEffectedEdge.begin() ; it != AllEffectedEdge.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        }
        
//----> change in global variables
        // volume and total area coupling
        double Volume_part_old = 0;
        double Area_part_old = 0;
        
        
        std::cout<<" here should be modefied 0092-3\n";

      /*  if((m_pState->GetVolumeCoupling())->GetState()|| (m_pState->GetApply_Constant_Area())->GetState()){
           
            if((m_pState->GetVolumeCoupling())->GetState()){
                for (std::vector<triangle *>::iterator it = Ctriangles.begin() ; it != Ctriangles.end(); ++it){
                    Volume_part_old+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(*it);
                    Area_part_old+=(*it)->GetArea();
                }
            }
            else{
                for (std::vector<triangle *>::iterator it = Ctriangles.begin() ; it != Ctriangles.end(); ++it){
                    Area_part_old+=(*it)->GetArea();
                }
            }
        }*/
        //- global curvature
        double vertex_area = 0;
        double vertex_Carea = 0;
        std::cout<<" threeedge should be fixed 1\n";
       /* if(m_pState->GetGlobalCurvature()->GetState()){
            double a11= ((it->PT1).pv1)->GetArea();
            double a12= ((it->PT1).pv2)->GetArea();
            double a13= ((it->PT1).pv3)->GetArea();
            double a21= ((it->PT2).pv1)->GetArea();
            double a22=((it->PT2).pv2)->GetArea();
            double a23= ((it->PT2).pv3)->GetArea();
            double c11=(((it->PT1).pv1)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c12=(((it->PT1).pv2)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c13=(((it->PT1).pv3)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c21=(((it->PT2).pv1)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c22=(((it->PT2).pv2)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c23=(((it->PT2).pv3)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            vertex_area=-(a11+a12+a13+a21+a22+a23);
            vertex_Carea=-(c11*a11+c12*a12+c13*a13+c21*a21+c22*a22+c23*a23);
        }*/
//----> end of: change in global variables should come in


            std::vector <triangle *> pair_tri = DoAScission(*it);

         
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv1);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv2);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv3);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv1);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv2);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv3);
            

//=== inclusion interaction energy;
            std::vector <links *> nl11 = ((it->PT1).pv1)->GetVLinkList();
            std::vector <links *> nl12 = ((it->PT1).pv2)->GetVLinkList();
            std::vector <links *> nl13 = ((it->PT1).pv3)->GetVLinkList();
            std::vector <links *> nl21 = ((it->PT2).pv1)->GetVLinkList();
            std::vector <links *> nl22 = ((it->PT2).pv2)->GetVLinkList();
            std::vector <links *> nl23 = ((it->PT2).pv3)->GetVLinkList();
            std::vector <links *> newAllEffectedEdge = nl11;
            newAllEffectedEdge.insert(newAllEffectedEdge.end(), nl12.begin(), nl12.end());
            newAllEffectedEdge.insert(newAllEffectedEdge.end(), nl13.begin(), nl13.end());
            newAllEffectedEdge.insert(newAllEffectedEdge.end(), nl21.begin(), nl21.end());
            newAllEffectedEdge.insert(newAllEffectedEdge.end(), nl22.begin(), nl22.end());
            newAllEffectedEdge.insert(newAllEffectedEdge.end(), nl23.begin(), nl23.end());
            RemoveFromLinkList((it->PT1).pl1,newAllEffectedEdge);
            RemoveFromLinkList((it->PT1).pl2,newAllEffectedEdge);
            RemoveFromLinkList((it->PT1).pl3,newAllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl1,newAllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl2,newAllEffectedEdge);
            RemoveFromLinkList((it->PT2).pl3,newAllEffectedEdge);
            
            for (std::vector<links *>::iterator it = newAllEffectedEdge.begin() ; it != newAllEffectedEdge.end(); ++it){
                enew+=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            }
//=== end of inclusion interaction energy

//----> change in global variables should come in
                double Volume_part_new = 0;
                double Area_part_new = 0;
        std::cout<<" here should be modefied 775683\n";

              /*  if((m_pState->GetVolumeCoupling())->GetState()|| (m_pState->GetApply_Constant_Area())->GetState()){
                   
                    if((m_pState->GetVolumeCoupling())->GetState()){
                        for (std::vector<triangle *>::iterator it = pair_tri.begin() ; it != pair_tri.end(); ++it){
                            Volume_part_new+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(*it);
                            Area_part_new+=(*it)->GetArea();
                        }
                    }
                    else{
                        for (std::vector<triangle *>::iterator it = pair_tri.begin() ; it != pair_tri.end(); ++it){
                            Area_part_new+=(*it)->GetArea();
                        }
                    }
                }*/
        
        //- global curvature
        double en_g_curve = 0;
        std::cout<<" threeedge should be fixed 2\n";

      /*  if(m_pState->GetGlobalCurvature()->GetState()){

            double a11= ((it->PT1).pv1)->GetArea();
            double a12= ((it->PT1).pv2)->GetArea();
            double a13= ((it->PT1).pv3)->GetArea();
            double a21= ((it->PT2).pv1)->GetArea();
            double a22=((it->PT2).pv2)->GetArea();
            double a23= ((it->PT2).pv3)->GetArea();
            double c11=(((it->PT1).pv1)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c12=(((it->PT1).pv2)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c13=(((it->PT1).pv3)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c21=(((it->PT2).pv1)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c22=(((it->PT2).pv2)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            double c23=(((it->PT2).pv3)->GetP1Curvature())+(((it->PT1).pv1)->GetP1Curvature());
            
            
            vertex_area+=(a11+a12+a13+a21+a22+a23);
            vertex_Carea+=(c11*a11+c12*a12+c13*a13+c21*a21+c22*a22+c23*a23);
            en_g_curve = m_pState->GetGlobalCurvature()->CalculateEnergyChange(vertex_area,vertex_Carea);

        }*/
            //--- change int energies of global
            double de_volume = 0;
            double de_global_area = 0;
        std::cout<<" here should be modefied 3333\n";

           /* if((m_pState->GetVolumeCoupling())->GetState()==true)
                de_volume = (m_pState->GetVolumeCoupling())->GetEnergyChange(step,Area_part_old,Volume_part_old,Area_part_new,Volume_part_new);
        
            if((m_pState->GetApply_Constant_Area())->GetState()==true)
                de_global_area =(m_pState->GetApply_Constant_Area())->GetEnergyChange(step,Area_part_old,Area_part_new);*/
//------ end of: change in global variable
  
            
            double de = enew - eold;
            double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
            if(de+de_volume+de_global_area+en_g_curve<0 || exp(-m_Beta*(de+de_volume+de_global_area+en_g_curve))>thermal){
                //--- Accepted
                m_pState->GetEnergyCalculator()->AddToTotalEnergy(de);
                
                std::cout<<" here should be modefied 1234\n";
                /*
                 
                (m_pState->GetVolumeCoupling())->UpdateArea_Volume(Area_part_old,Volume_part_old,Area_part_new,Volume_part_new);
                (m_pState->GetApply_Constant_Area())->UpdateArea(Area_part_old,Area_part_new);
                //-- global curvature update of energy
                if(m_pState->GetGlobalCurvature()->GetState()){
                    m_pState->GetGlobalCurvature()->UpdateEnergyChange(vertex_area,vertex_Carea);

                }*/
            }
            else{ // reject the move
                ReverseAScission(*it, pair_tri[0], pair_tri[1]);
                double ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv1);
                ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv2);
                ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv3);
                ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv1);
                ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv2);
                ten= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv3);
                //--- only the new is enough as the old (coonecting edges are not affected, they were ghost)
                for (std::vector<links *>::iterator it = newAllEffectedEdge.begin() ; it != newAllEffectedEdge.end(); ++it){
                    enew+=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                }
            }
    }
   // auto end = std::chrono::steady_clock::now();
   // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
   // std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
    return true;
}
// it creates a triangle and place it to the active trinagles list
triangle * Three_Edge_Scission::CreateATriangleFromAPotentialTriangle(pot_triangle &p1)
{
//---> initate the copy of the three links for the potential triangle
    (p1.pl1)->SetCopy();
    (p1.pl2)->SetCopy();
    (p1.pl3)->SetCopy();
    
//---> create the triangle
    //--- select the last ghost triangle
    triangle *gt1 = m_pGhostT[m_pGhostT.size()-1];
    //--- update its vertex
    gt1->UpdateVertex(p1.pv1,p1.pv2,p1.pv3);

    //--- put the triangle from ghost to active
    ( m_pActiveT).push_back(gt1);
    m_pGhostT.pop_back();
    
//---> now make the changes in the edges
    //--- update their triangles
    (p1.pl1)->UpdateTriangle(gt1);
    (p1.pl2)->UpdateTriangle(gt1);
    (p1.pl3)->UpdateTriangle(gt1);
    //---- update their next two edges
    p1.pl1->UpdateNeighborLink1(p1.pl2);
    p1.pl1->UpdateNeighborLink2(p1.pl3);
    p1.pl2->UpdateNeighborLink1(p1.pl3);
    p1.pl2->UpdateNeighborLink2(p1.pl1);
    p1.pl3->UpdateNeighborLink1(p1.pl1);
    p1.pl3->UpdateNeighborLink2(p1.pl2);
    //--- update their v3 (v2 and v1 remain the same)
    p1.pl1->UpdateV3(p1.pv3);
    p1.pl2->UpdateV3(p1.pv1);
    p1.pl3->UpdateV3(p1.pv2);
    
//---> update the vertices triangle list
    p1.pv1->AddtoTraingleList(gt1);
    p1.pv2->AddtoTraingleList(gt1);
    p1.pv3->AddtoTraingleList(gt1);
    
//---> update geometry; we cannot update the geometry of the the vertices because they have some edges that should be removed
    //--- update the normal and area of the triangle
    gt1->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be obtained
    //--- update the three links notmal and shape Operator
    (p1.pl1)->UpdateNormal();
    (p1.pl2)->UpdateNormal();
    (p1.pl3)->UpdateNormal();
    (p1.pl1)->UpdateShapeOperator(m_pBox);
    (p1.pl2)->UpdateShapeOperator(m_pBox);
    (p1.pl3)->UpdateShapeOperator(m_pBox);

    return gt1;
}
// this function cuts the neck made of p1 and p2, i.e., pair
std::vector <triangle *> Three_Edge_Scission::DoAScission(pair_pot_triangle &pair){
    
    std::vector <triangle *> createdtriangles;
    if(m_pGhostT.size()<2){
        std::cout<<" ---> [not an error] not enough reserved trinagles; you may restart the simulation "<<std::endl;
        exit(0);
    }
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
//---> generate the two triangles and store them in the return vector
    //--- create the triangles
    createdtriangles.push_back(CreateATriangleFromAPotentialTriangle(p1));
    createdtriangles.push_back(CreateATriangleFromAPotentialTriangle(p2));

//---> remove the links that connect the two potential trinagles; we just remove them from different lists and send them to ghost
    //-- these edges are stored in the Clinks that are obtained in when the pair is created
    //-- note, the mirror links do not exist in this list, also their triangles should be deleted
    //-- these two triangles now have been created
    std::vector <links *> Clinks = pair.ConnectingLinks;       // this does not include the mirror links
    for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){
        //-- remove the edge and its mirror
        RemoveFromLinkList(*it,  m_pActiveL);
        RemoveFromLinkList(*it,  m_pHL);
        RemoveFromLinkList(*it,  m_pMHL);
        RemoveFromLinkList((*it)->GetMirrorLink(),  m_pActiveL);
        RemoveFromLinkList((*it)->GetMirrorLink(),  m_pHL);
        RemoveFromLinkList((*it)->GetMirrorLink(),  m_pMHL);
        AddtoVectorCarefully(*it,m_pGhostL);
        AddtoVectorCarefully((*it)->GetMirrorLink(),m_pGhostL);
        //-- remove the edge from the list of the vertices
        ((*it)->GetV1())->RemoveFromLinkList(*it);
        ((*it)->GetV2())->RemoveFromLinkList((*it)->GetMirrorLink());
        //-- remove the assosciated triangle
        RemoveFromTriangleList((*it)->GetTriangle(), m_pActiveT);
        AddtoVectorCarefully((*it)->GetTriangle(),m_pGhostT);
        //-- remove the associated triangle from the vertices list
        ((*it)->GetTriangle()->GetV1())->RemoveFromTraingleList((*it)->GetTriangle());
        ((*it)->GetTriangle()->GetV2())->RemoveFromTraingleList((*it)->GetTriangle());
        ((*it)->GetTriangle()->GetV3())->RemoveFromTraingleList((*it)->GetTriangle());
        //-- make sure that  v1 and v2 are not connected by v-list
        ((*it)->GetV1())->RemoveFromNeighbourVertex((*it)->GetV2());
        ((*it)->GetV2())->RemoveFromNeighbourVertex((*it)->GetV1());
    } // end of     "for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){"

//--> update geometry
    //-- update geometry of the 6 vertices
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv3);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv3);

    return createdtriangles;
}
// this is the exact reverse action of DoAScission; different from DoAFussion
bool Three_Edge_Scission::ReverseAScission(pair_pot_triangle &pair , triangle *t1, triangle *t2){

//--- getting p1 and p2
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
//-- we take the triangle3 t1 and t2 to the ghost and remove them from associated vertices and edges
    //-- send t1 and t2 to the ghost
    RemoveFromTriangleList(t1, m_pActiveT);
    RemoveFromTriangleList(t2, m_pActiveT);
    m_pGhostT.push_back(t1);
    m_pGhostT.push_back(t2);
    //--- remove t1 and t2 from their v->t list
    (t1->GetV1())->RemoveFromTraingleList(t1);
    (t1->GetV2())->RemoveFromTraingleList(t1);
    (t1->GetV3())->RemoveFromTraingleList(t1);
    (t2->GetV1())->RemoveFromTraingleList(t2);
    (t2->GetV2())->RemoveFromTraingleList(t2);
    (t2->GetV3())->RemoveFromTraingleList(t2);
    //--- reverse the edges to previous value
    (p1.pl1)->Reverse2PreviousCopy();
    (p1.pl2)->Reverse2PreviousCopy();
    (p1.pl3)->Reverse2PreviousCopy();
    (p2.pl1)->Reverse2PreviousCopy();
    (p2.pl2)->Reverse2PreviousCopy();
    (p2.pl3)->Reverse2PreviousCopy();

     //
//==== from the pair, we have all the links that were used to connect t1 and t2, we can recomver it well
    std::vector <links *> Clinks = pair.ConnectingLinks;       // this does not include the mirror links
    for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){

        //--removing it and its mirror from ghost and add to real containors
        AddtoVectorCarefully((*it), m_pActiveL);
        AddtoVectorCarefully((*it), m_pMHL);
        AddtoVectorCarefully((*it)->GetMirrorLink(), m_pActiveL);
        AddtoVectorCarefully((*it)->GetMirrorLink(), m_pHL);
        RemoveFromLinkList(*it,m_pGhostL);
        RemoveFromLinkList((*it)->GetMirrorLink(),m_pGhostL);
        //-- also bring the triangles to live again
        RemoveFromTriangleList((*it)->GetTriangle(),m_pGhostT);
        AddtoVectorCarefully((*it)->GetTriangle(), m_pActiveT);
        //-- add the links to the vertex linklist
        ((*it)->GetV1())->AddtoLinkListCarefully((*it));
        ((*it)->GetV2())->AddtoLinkListCarefully((*it)->GetMirrorLink());
        //--- add the triangles to the vertex list
        ((*it)->GetV1())->AddtoTriangleListCarefully((*it)->GetTriangle());
        ((*it)->GetV2())->AddtoTriangleListCarefully((*it)->GetTriangle());
        ((*it)->GetV3())->AddtoTriangleListCarefully((*it)->GetTriangle());
        //--- make v1 and v2 nighbour again
        ((*it)->GetV1())->AddtoNeighbourVertexCarefully((*it)->GetV2());
        ((*it)->GetV2())->AddtoNeighbourVertexCarefully((*it)->GetV1());
        //--- note l1 and its mirror were as their old as never changed
    } //   End of "for (std::vector<links*>::iterator it = Clinks.begin() ; it != Clinks.end(); it++){"
//--> update geometry
    //-- update geometry of the 6 vertices
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p1.pv3);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv1);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv2);
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(p2.pv3);
    
    return true;
}
//== this function get a mesh and search through it and finds possible fission sites
std::vector<pair_pot_triangle> Three_Edge_Scission::FindPotentialTriangles(){

//---- this part searches through all the links and finds possible triangles that do not exist. This means, it finds triple vertices (v1,v2,v3) that are connected by edges but such
// trinagle is not defined.
    
    int id = 0;
    std::vector<pot_triangle> list;
    std::vector<links*>  all_link =  m_pActiveL;
    for (std::vector<links*>::iterator il1 = all_link.begin() ; il1 != all_link.end(); il1++){
        vertex *pv1 = (*il1)->GetV1();
        vertex *pv2 = (*il1)->GetV2();
        std::vector<links*>  v2_l = pv2->GetVLinkList();
        for (std::vector<links*>::iterator il2 = v2_l.begin() ; il2 != v2_l.end(); il2++){
            vertex *pv3 = (*il2)->GetV2();
            std::vector<links*>  v3_l = pv3->GetVLinkList();
                for (std::vector<links*>::iterator il3 = v3_l.begin() ; il3 != v3_l.end(); il3++){
                    if(pv1==(*il3)->GetV2() && (*il3)->GetV3()!=pv2 && ((*il3)->GetMirrorLink())->GetV3()!=pv2 ){
                        if( pv2->m_VertexType==0 && pv3->m_VertexType==0 && pv1->m_VertexType==0){
                            // pv1, pv2, pv3 are our vertices
                                                    
                            pot_triangle PotT;
                            PotT.id = id; id++;
                            PotT.cid = -1;
                            PotT.pv1 = pv1; PotT.pv2 = pv2; PotT.pv3 = pv3;
                            PotT.pl1= (*il1); PotT.pl2= (*il2); PotT.pl3= (*il3);
                            list.push_back(PotT);
                            //--- we remove these links from all link containor so we do not search through them again
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il1)->GetMirrorLink()), all_link.end());

                        }
                    }
            }
        }
    }
//===================
    std::vector <pair_pot_triangle> pair_list;
    int idpair=0;
for (int i=0;i<list.size();i++)
{
    for (int j=i+1;j<list.size();j++){
        pair_pot_triangle temp = connected_2pot_triangles(list[i],list[j]);
        if(temp.state){
            temp.id = idpair;
            idpair++;
            pair_list.push_back(temp);
        }
    }
}
    
    return pair_list;
}
//-- connected_2pot_triangles function check if the two potential trinagles are well connected for a fission. not, T1 and T2 do not exist yet, but they can apear if v1,v2,v3 get
//          v1------------v4            disconnected from v4, v5, v6. This function checks for such cases
//         /T1\         / T2\
//        v2--v3-------v4---v6
pair_pot_triangle Three_Edge_Scission::connected_2pot_triangles(pot_triangle potT1, pot_triangle potT2)
{
    pair_pot_triangle ReturningPair;
    ReturningPair.state = false;
    
    if(potT2.cid != -1 || potT1.cid != -1)
        return ReturningPair;
    if(potT2.pv1 == potT1.pv1 || potT2.pv1 == potT1.pv2 || potT2.pv1 == potT1.pv3)  // if the two trinagles have a shared vertex
        return ReturningPair;
    if(potT2.pv2 == potT1.pv1 || potT2.pv2 == potT1.pv2 || potT2.pv2 == potT1.pv3) // if the two trinagles have a shared vertex
        return ReturningPair;
    if(potT2.pv3 == potT1.pv1 || potT2.pv3 == potT1.pv2 || potT2.pv3 == potT1.pv3) // if the two trinagles have a shared vertex
        return ReturningPair;
    
    int test1 = 0;
    int test2 = 0;
    int test3 = 0;
    std::vector <vertex *> nv1 = (potT1.pv1)->GetVNeighbourVertex();
    std::vector <vertex *> nv2 = (potT1.pv2)->GetVNeighbourVertex();
    std::vector <vertex *> nv3 = (potT1.pv3)->GetVNeighbourVertex();
    int connect[3][3] = {0};   // a matrix that shows how each vertex of the two trinagle are connected
    // here we attempt to fill the connect matrix and make sure that for each vertex from T1 at least one connected vertex from T2 exist
    // this condidtion is stored in test1-3 varible. However, it can happen in which all the verices from T1 is connected to one vertex from T2
    // therefore we check the connect matrix
    for (std::vector<vertex*>::iterator it = nv1.begin() ; it != nv1.end(); it++){
        
        if(potT2.pv1 == (*it) || potT2.pv2 == (*it) || potT2.pv3 == (*it))
            test1 = true;
        if(potT2.pv1 == (*it))
            connect[0][0]=1;
        if(potT2.pv2 == (*it))
            connect[1][0]=1;
        if(potT2.pv3 == (*it))
            connect[2][0]=1;
        
    }
    for (std::vector<vertex*>::iterator it = nv2.begin() ; it != nv2.end(); it++){
       
        if(potT2.pv1 == (*it) || potT2.pv2 == (*it) || potT2.pv3 == (*it))
            test2 = true;
        if(potT2.pv1 == (*it))
            connect[0][1]=1;
        if(potT2.pv2 == (*it))
            connect[1][1]=1;
        if(potT2.pv3 == (*it))
            connect[2][1]=1;
    }
    for (std::vector<vertex*>::iterator it = nv3.begin() ; it != nv3.end(); it++){
        
        if(potT2.pv1 == (*it) || potT2.pv2 == (*it) || potT2.pv3 == (*it))
            test3 = true;
        if(potT2.pv1 == (*it))
            connect[0][2]=1;
        if(potT2.pv2 == (*it))
            connect[1][2]=1;
        if(potT2.pv3 == (*it))
            connect[2][2]=1;
    }

//-- only if the trinagles are connected; we check if T2 is also well connected to T1. connect mat is to avoid repeataion of previous process
    if(test1==true && test2==true && test3==true )
    {
        for (int i=0;i<3;i++){
            int row = 0;
            for (int j=0;j<3;j++){
                row+= connect[i][j];
            }
            if(row==0)
                return ReturningPair;
        }
        
        potT2.cid = potT1.id;
        potT1.cid = potT2.id;
    }
    else{
        return ReturningPair;
    }
    
    //--- here we should find all the links that connects the  vertices of P1 to P2.
    // how to return this
    if(CorrectOrientation(potT1,potT2) && CorrectOrientation(potT2,potT1)){ // this means that the p1 edges should be so that the removal edge disconect
                                                                 //vertices of p2 from p1, not any other vertices
        ReturningPair.PT1 = potT1;
        ReturningPair.PT2 = potT2;
        std::vector <links *> Clinks;
        std::vector <triangle *> Ctriangles;
        ReturningPair.state = true;

        // Iterate through the links connected to the first vertex of potT1
        std::vector<links*> n_p1links = (potT1.pv1)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p1links.begin(); it != n_p1links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.push_back((*it)->GetTriangle());
            }
        }

        // Iterate through the links connected to the second vertex of potT1
        std::vector<links*> n_p2links = (potT1.pv2)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p2links.begin(); it != n_p2links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.push_back((*it)->GetTriangle());
            }
        }

        // Iterate through the links connected to the third vertex of potT1
        std::vector<links*> n_p3links = (potT1.pv3)->GetVLinkList();
        for (std::vector<links*>::iterator it = n_p3links.begin(); it != n_p3links.end(); ++it) {
            if ((*it)->GetV2() == potT2.pv1 || (*it)->GetV2() == potT2.pv2 || (*it)->GetV2() == potT2.pv3) {
                Clinks.push_back(*it);
                Ctriangles.push_back((*it)->GetTriangle());
            }
        }
        ReturningPair.ConnectingLinks = Clinks;
        ReturningPair.ConnectingTriangles = Ctriangles;
    }
    else
    {
        std::cout<<"---> error: this should not happen \n";
    }
    
    
    
    return ReturningPair;
}
//--- v1,v2 and v3 may not have the correct trinagluation Orientation.
 //      >v         ;correct orientation means to disconnect links connected to vertices in the pot_triangle 2
//  l1  /  \>  l2   ;because for each triple of v1, v2 and v3, there is two ways to create a trinagle but each will
//    v1<---v3      ; leads to removal of diffierent edges.
bool Three_Edge_Scission::CorrectOrientation(pot_triangle &p1,pot_triangle &p2)
{
    links* ml1 = (p1.pl1)->GetMirrorLink();
    if((p1.pl1)->GetV3()==p2.pv1 || (p1.pl1)->GetV3()==p2.pv2 || (p1.pl1)->GetV3()==p2.pv3){
        // Orientation is correct, no change is needed
        return true;
    }
    else if(ml1->GetV3()==p2.pv1 || ml1->GetV3()==p2.pv2 || ml1->GetV3()==p2.pv3){
        // Orientation is not correct, we reverse it
        links* ml2 = (p1.pl2)->GetMirrorLink();
        links* ml3 = (p1.pl3)->GetMirrorLink();
        vertex *v1 = p1.pv1;
        vertex *v2 = p1.pv2;
        p1.pv1 = v2;
        p1.pv2 = v1;
        
        p1.pl1 = ml1;
        p1.pl2 = ml3;
        p1.pl3 = ml2;


        return true;
    }
    else {
        return false;
    }
    
    return true;
}

bool Three_Edge_Scission::DoAFussion(pair_pot_triangle pair)
{
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
    
    return false;
}
void Three_Edge_Scission::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void Three_Edge_Scission::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
template<typename T>
bool Three_Edge_Scission::AddtoVectorCarefully(T* item, std::vector<T*>& vect) {
    // Check if the item already exists in the list
    for (typename std::vector<T*>::iterator it = vect.begin(); it != vect.end(); ++it) {
        if (*it == item)
            return false;
    }
    // Item does not exist, add it to the list
    vect.push_back(item);
    return true;
}
// this should be deleted at the end

template<typename T>
void Three_Edge_Scission::KeepOneOccurrence(std::vector<T*> &vec){
   
    std::sort(vec.begin(), vec.end()); // Sort the vector
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); // Remove duplicates
}
std::string Three_Edge_Scission::CurrentState(){
    
    std::string state = AbstractDynamicTopology::GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}



