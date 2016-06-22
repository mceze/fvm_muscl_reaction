/*
 *  Init.c
 *  fvm
 *
 *  Created by Marco Ceze on 8/13/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "tools_IO.h"
#include "tools_memory.h"
#include "structure_all.h"
#include "structure_data.h"
#include "structure_IO.h"
#include "tools_solver.h"

/****************************************************************************/
//Function fvm_Init
int fvm_Init(fvm_All *All)
{
  int i, j, k, l, d, egrp, elem, node, fgrp, face;
  int nodeA, nodeB, nodeC, nodeD;
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  fvm_IO *IO;
  
  Data = All->Data;
  Mesh = All->Mesh;
  IO = All->IO;
  
  Data->CurrentTime = 0.0;
  
  /****************************************************************************/
  //allocating the arrays that store the data
  
  //Element Group related data
  fvm_Alloc((void **)&Data->ElemGroup, Mesh->nElemGroup, sizeof(fvm_EGrpData));
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    fvm_Alloc2((void ***)&Data->ElemGroup[egrp].ElemCenter, 
               Mesh->ElemGroup[egrp].nElem, Mesh->Dim, sizeof(double));
    fvm_Alloc((void **)&Data->ElemGroup[egrp].Dt, 
              Mesh->ElemGroup[egrp].nElem, sizeof(double));
    fvm_Alloc((void **)&Data->ElemGroup[egrp].State, 
              Mesh->ElemGroup[egrp].nElem, sizeof(fvm_State));
    fvm_Alloc2((void ***)&Data->ElemGroup[egrp].dU_dt, 
               Mesh->ElemGroup[egrp].nElem, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
    fvm_Alloc((void **)&Data->ElemGroup[egrp].Volume, 
              Mesh->ElemGroup[egrp].nElem, sizeof(double));
    if (IO->FluxType == JAMESON){
      fvm_Alloc2((void ***)&Data->ElemGroup[egrp].C, 
                Mesh->ElemGroup[egrp].nElem, Mesh->Dim+2, sizeof(double));
      fvm_Alloc2((void ***)&Data->ElemGroup[egrp].D, 
                Mesh->ElemGroup[egrp].nElem, Mesh->Dim+2, sizeof(double));
      fvm_Alloc2((void ***)&Data->ElemGroup[egrp].Nabla2_Q, 
                 Mesh->ElemGroup[egrp].nElem, Mesh->Dim+2, sizeof(double));
      fvm_Alloc((void **)&Data->ElemGroup[egrp].A, 
                Mesh->ElemGroup[egrp].nElem, sizeof(double));
      fvm_Alloc((void **)&Data->ElemGroup[egrp].nu, 
                Mesh->ElemGroup[egrp].nElem, sizeof(double));
    }
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      fvm_Alloc((void **)&Data->ElemGroup[egrp].State[elem].U, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
      fvm_Alloc((void **)&Data->ElemGroup[egrp].State[elem].W, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
    }
  }
  //Boundary Group related data
  fvm_Alloc((void **)&Data->BFaceGroup, Mesh->nBFaceGroup, sizeof(fvm_BFGrpData));
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    fvm_Alloc((void **)&Data->BFaceGroup[fgrp].Area, 
              Mesh->BFaceGroup[fgrp].nFace, sizeof(double));
    fvm_Alloc2((void ***)&Data->BFaceGroup[fgrp].Normal, 
               Mesh->BFaceGroup[fgrp].nFace, Mesh->Dim, sizeof(double));
    fvm_Alloc2((void ***)&Data->BFaceGroup[fgrp].Tangent, 
               Mesh->BFaceGroup[fgrp].nFace, Mesh->Dim, sizeof(double));
    fvm_Alloc((void **)&Data->BFaceGroup[fgrp].State, 
              Mesh->BFaceGroup[fgrp].nFace, sizeof(fvm_State));
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
      fvm_Alloc((void **)&Data->BFaceGroup[fgrp].State[face].W, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
      fvm_Alloc((void **)&Data->BFaceGroup[fgrp].State[face].U, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
    }
  }
  //Interior faces related data
  fvm_Alloc((void **)&Data->IFace.Area, 
            Mesh->nIFace, sizeof(double));
  fvm_Alloc2((void ***)&Data->IFace.Normal, 
             Mesh->nIFace, Mesh->Dim, sizeof(double));
  fvm_Alloc2((void ***)&Data->IFace.Tangent, 
             Mesh->nIFace, Mesh->Dim, sizeof(double));
  
  /****************************************************************************/
  //Initializing State and Volume on Elements
  for (j = 0; j < IO->nInitStates; j++){
    fvm_ConvertState(IO->InitState[j], IO->InitType[j], Mesh->Dim, IO->nChemSpecies);
    for (i = 0; i < IO->nEGroups[j]; i++){
      egrp = IO->InitGroups[j][i];
      for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
          Data->ElemGroup[egrp].State[elem].U[k] = IO->InitState[j].U[k];
          Data->ElemGroup[egrp].State[elem].W[k] = IO->InitState[j].W[k];
        }
        for (d = 0; d < Mesh->Dim; d++){
          Data->ElemGroup[egrp].ElemCenter[elem][d] = 0.0;
          for (l = 0; l < Mesh->ElemGroup[egrp].nNode; l++){
            node = Mesh->ElemGroup[egrp].Node[elem][l];
            Data->ElemGroup[egrp].ElemCenter[elem][d] += Mesh->Coord[node][d];
          }
          Data->ElemGroup[egrp].ElemCenter[elem][d] = Data->ElemGroup[egrp].ElemCenter[elem][d]/Mesh->ElemGroup[egrp].nNode;
        }
        Data->ElemGroup[egrp].Dt[elem] = 0.0;
        nodeA = Mesh->ElemGroup[egrp].Node[elem][0];
        nodeB = Mesh->ElemGroup[egrp].Node[elem][1];
        nodeC = Mesh->ElemGroup[egrp].Node[elem][2];
        nodeD = Mesh->ElemGroup[egrp].Node[elem][3];
        
        Data->ElemGroup[egrp].Volume[elem] = 0.5*((Mesh->Coord[nodeB][0]-Mesh->Coord[nodeD][0])
                                                  *(Mesh->Coord[nodeC][1]-Mesh->Coord[nodeA][1])
                                                  -(Mesh->Coord[nodeC][0]-Mesh->Coord[nodeA][0])
                                                  *(Mesh->Coord[nodeB][1]-Mesh->Coord[nodeD][1]));
      }
    }
  }
  
  //Initializing Face data
  fvm_CalculateFaceData(All);
  
  return 0;
}

/****************************************************************************/
//Function fvm_CalculateFaceData
int fvm_CalculateFaceData(fvm_All *All)
{
  int fgrp, face, nodes[2];//only works for 2D
  int faceL, elemL, egrpL;
  double temp;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  
  Mesh = All->Mesh;
  Data = All->Data;
  
  //boundary faces
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
      egrpL = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
      elemL = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
      faceL = Mesh->BFaceGroup[fgrp].BFace[face].Face;
      if (faceL == Mesh->ElemGroup[egrpL].nFace - 1){
        nodes[0] = Mesh->ElemGroup[egrpL].Node[elemL][Mesh->ElemGroup[egrpL].nNode-1];
        nodes[1] = Mesh->ElemGroup[egrpL].Node[elemL][0];
      }
      else{
        nodes[0] = Mesh->ElemGroup[egrpL].Node[elemL][faceL];
        nodes[1] = Mesh->ElemGroup[egrpL].Node[elemL][faceL + 1];
      }
      temp = pow(Mesh->Coord[nodes[0]][0]-Mesh->Coord[nodes[1]][0],
                 2.0)+pow(Mesh->Coord[nodes[0]][1]-Mesh->Coord[nodes[1]][1],2.0);
      Data->BFaceGroup[fgrp].Area[face] = sqrt(temp);
      
      Data->BFaceGroup[fgrp].Normal[face][0] = (Mesh->Coord[nodes[1]][1]-Mesh->Coord[nodes[0]][1])/Data->BFaceGroup[fgrp].Area[face];
      Data->BFaceGroup[fgrp].Normal[face][1] = -(Mesh->Coord[nodes[1]][0]-Mesh->Coord[nodes[0]][0])/Data->BFaceGroup[fgrp].Area[face];
      Data->BFaceGroup[fgrp].Tangent[face][0] = -Data->BFaceGroup[fgrp].Normal[face][1];
      Data->BFaceGroup[fgrp].Tangent[face][1] = Data->BFaceGroup[fgrp].Normal[face][0];
    }
  }
  
  //internal faces
  
  for (face = 0; face < Mesh->nIFace; face++){
    egrpL = Mesh->IFace[face].ElemGroupL;
    elemL = Mesh->IFace[face].ElemL;
    faceL = Mesh->IFace[face].FaceL;
    if (faceL == Mesh->ElemGroup[egrpL].nFace - 1){
      nodes[0] = Mesh->ElemGroup[egrpL].Node[elemL][Mesh->ElemGroup[egrpL].nNode-1];
      nodes[1] = Mesh->ElemGroup[egrpL].Node[elemL][0];
    }
    else{
      nodes[0] = Mesh->ElemGroup[egrpL].Node[elemL][faceL];
      nodes[1] = Mesh->ElemGroup[egrpL].Node[elemL][faceL + 1];
    }
    temp = pow(Mesh->Coord[nodes[0]][0]-Mesh->Coord[nodes[1]][0],
               2.0)+pow(Mesh->Coord[nodes[0]][1]-Mesh->Coord[nodes[1]][1],2.0);
    Data->IFace.Area[face] = sqrt(temp);
    Data->IFace.Normal[face][0] = (Mesh->Coord[nodes[1]][1]-Mesh->Coord[nodes[0]][1])/Data->IFace.Area[face];
    Data->IFace.Normal[face][1] = -(Mesh->Coord[nodes[1]][0]-Mesh->Coord[nodes[0]][0])/Data->IFace.Area[face];
    Data->IFace.Tangent[face][0] = -Data->IFace.Normal[face][1];
    Data->IFace.Tangent[face][1] = Data->IFace.Normal[face][0];
  }
  
  return 0;
}



