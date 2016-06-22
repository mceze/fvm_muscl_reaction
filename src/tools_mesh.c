/*
 *  tools_mesh.c
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
//Function fvm_NeighborAcrossFace
int fvm_NeighborAcrossFace(fvm_Mesh *Mesh, int egrp, int elem, int face, int *neighbor_grp, int *neighbor, int *neighbor_face, int *is_it_boundary)
{
  int f;
  //check if it a boundary
  if (Mesh->ElemGroup[egrp].Face[elem][face].Group != INTERIORFACE){
    *is_it_boundary = 1;
    return 0;
  }
  else
    *is_it_boundary = 0;
  
  //it is an interior face
  f = Mesh->ElemGroup[egrp].Face[elem][face].Number;
  if (Mesh->IFace[f].ElemGroupL != egrp || Mesh->IFace[f].ElemL != elem){
    *neighbor_grp = Mesh->IFace[f].ElemGroupL;
    *neighbor = Mesh->IFace[f].ElemL;
    *neighbor_face = Mesh->IFace[f].FaceL;
  }
  else {
    *neighbor_grp = Mesh->IFace[f].ElemGroupR;
    *neighbor = Mesh->IFace[f].ElemR;
    *neighbor_face = Mesh->IFace[f].FaceR;
  }
  
  return 0;
}

/****************************************************************************/
//Function fvm_FindOppositeFace
int fvm_FindOppositeFace(int face, int egrp, fvm_All *All)
{
  int OppositeFace, nFace;
  nFace = All->Mesh->ElemGroup[egrp].nFace;
  
  //for now it works only for quads
  OppositeFace = face + 2;
  if (OppositeFace >= nFace)
    OppositeFace = OppositeFace - nFace;
  
  return OppositeFace;
}
