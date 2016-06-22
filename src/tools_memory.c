/*
 *  tools_memory.c
 *  fvm
 *
 *  Created by Marco Ceze on 7/22/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "tools_memory.h"
#include "structure_mesh.h"
#include "structure_data.h"
#include "structure_IO.h"
#include "structure_all.h"

/****************************************************************************/
//Function Alloc
int fvm_Alloc(void **pchunk, int n, int size)
{
  int totalsize;
  
  (*pchunk) = NULL;
  
  totalsize = n*size;
  
  (*pchunk) = (void *)malloc(totalsize);
  
  return 0;
}

/****************************************************************************/
//Function Alloc2
int fvm_Alloc2( void ***pchunk, int n1, int n2, int size)
{
  char *temp;
  int i;
  int totalsize;
  
  (*pchunk) = NULL;
  
  totalsize = n1*n2*size;
  
  temp = (char *)malloc( totalsize);
  
  (*pchunk) = (void **)malloc( n1*sizeof(char *) );
  
  for(i = 0; i<n1; i++)
    (*pchunk)[i] = temp + i*n2*size;
  
  return 0;
}


/****************************************************************************/
//Function CreateMesh
int fvm_CreateMesh(fvm_Mesh **Mesh)
{    
  fvm_Alloc((void **)Mesh, 1, sizeof(fvm_Mesh));
  
  (*Mesh)->Dim   = 0;
  (*Mesh)->nNode = 0;
  (*Mesh)->Coord = NULL;
  (*Mesh)->nIFace = 0;
  (*Mesh)->IFace  = NULL;
  (*Mesh)->nBFaceGroup = 0;  
  (*Mesh)->BFaceGroup  = NULL;
  (*Mesh)->nElemGroup = 0;
  (*Mesh)->ElemGroup  = NULL;
  
  return 0;
}

/****************************************************************************/
//Function CreateData
int fvm_CreateData(fvm_Data **Data)
{    
  fvm_Alloc((void **)Data, 1, sizeof(fvm_Data));
  
  (*Data)->ElemGroup = NULL;
  (*Data)->BFaceGroup = NULL;
  
  return 0;
}

/****************************************************************************/
//Function CreateIO
int fvm_CreateIO(fvm_IO **IO)
{    
  fvm_Alloc((void **)IO, 1, sizeof(fvm_IO));
  
  (*IO)->nInitStates = 0;
  (*IO)->InitState = NULL;
  (*IO)->InitGroups = NULL;
  (*IO)->nEGroups = NULL;
  (*IO)->InitType = NULL;
  sprintf((*IO)->MeshFile, "mesh.neu");
  (*IO)->CFL = 0.0;
  (*IO)->itmax = 0;
  (*IO)->ConstDt = 1;  
  
  return 0;
}

/****************************************************************************/
//Function CreateAll
int fvm_CreateAll(fvm_All **pAll)
{
  fvm_All *All;
  
  fvm_Alloc((void **)pAll, 1, sizeof(fvm_All));
  All = (*pAll);
  
  fvm_CreateMesh(&(All->Mesh));
  fvm_CreateData(&(All->Data));
  fvm_CreateIO(&(All->IO));
  
  return 0;
}

/****************************************************************************/
//Function Release
void fvm_Release(void *chunk)
{
  if (chunk == NULL) return;
  free((void *)chunk);
}

/****************************************************************************/
//Function Release2
void fvm_Release2(void **chunk)
{
  if (chunk == NULL) return;
  free((void * ) chunk[0]);
  free((void **) chunk   );
}

/****************************************************************************/
//Function Release3
void fvm_Release3(void ***chunk)
{
  if (chunk == NULL) return;
  free( (void * )  chunk[0][0]);
  free( (void **)  chunk[0]   );
  free( (void ***) chunk   );
}

/****************************************************************************/
//Function DestroyMesh
int fvm_DestroyMesh(fvm_Mesh *Mesh)
{
  int i;
  if (Mesh == NULL) return 0;
  
  fvm_Release2((void **)Mesh->Coord);
  fvm_Release((void *)Mesh->IFace);
  for (i = 0; i < Mesh->nBFaceGroup; i++)
    fvm_Release((void *)Mesh->BFaceGroup[i].BFace);
  fvm_Release((void *)Mesh->BFaceGroup);
  for (i = 0; i < Mesh->nElemGroup; i++){
    fvm_Release2((void **)Mesh->ElemGroup[i].Node);
    fvm_Release2((void **)Mesh->ElemGroup[i].Face);
  }
  fvm_Release((void *)Mesh->ElemGroup);
  
  fvm_Release((void *)Mesh);
  
  return 0;
}

/****************************************************************************/
//Function DestroyData
int fvm_DestroyData(fvm_All *All)
{
  int i, j;
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  Data = All->Data;
  Mesh = All->Mesh;
  
  if (Data == NULL) return 0;
  for (i = 0; i < Mesh->nElemGroup; i++){
    fvm_Release2((void **)Data->ElemGroup[i].ElemCenter);
    fvm_Release((void *)Data->ElemGroup[i].Dt);
    for (j = 0; j < Mesh->ElemGroup[i].nElem; j++){
      fvm_Release((void *)Data->ElemGroup[i].State[j].W);
      fvm_Release((void *)Data->ElemGroup[i].State[j].U);
    }
    fvm_Release((void *)Data->ElemGroup[i].State);
    fvm_Release2((void **)Data->ElemGroup[i].dU_dt);
    fvm_Release((void *)Data->ElemGroup[i].Volume);
    if (All->IO->FluxType == JAMESON){
      fvm_Release2((void **)Data->ElemGroup[i].C);
      fvm_Release2((void **)Data->ElemGroup[i].D);
      fvm_Release2((void **)Data->ElemGroup[i].Nabla2_Q);
      fvm_Release((void *)Data->ElemGroup[i].A);
      fvm_Release((void *)Data->ElemGroup[i].nu);
    }
  }
  fvm_Release((void *)Data->ElemGroup);
  
  for (i = 0; i < Mesh->nBFaceGroup; i++){
    for (j = 0; j < Mesh->BFaceGroup[i].nFace; j++){
      fvm_Release((void *)Data->BFaceGroup[i].State[j].W);
      fvm_Release((void *)Data->BFaceGroup[i].State[j].U);
    }
    fvm_Release((void *)Data->BFaceGroup[i].State);
    fvm_Release((void *)Data->BFaceGroup[i].Area);
    fvm_Release2((void **)Data->BFaceGroup[i].Normal);
    fvm_Release2((void **)Data->BFaceGroup[i].Tangent);
  }
  fvm_Release((void *)Data->BFaceGroup);
  
  fvm_Release((void *)Data->IFace.Area);
  fvm_Release2((void **)Data->IFace.Normal);
  fvm_Release2((void **)Data->IFace.Tangent);
  
  fvm_Release((void *)Data);
  
  return 0;
}

/****************************************************************************/
//Function DestroyIO
int fvm_DestroyIO(fvm_IO *IO)
{
  int i;
  
  for (i = 0; i < IO->nInitStates; i++){
    fvm_Release((void *)IO->InitState[i].W);
    fvm_Release((void *)IO->InitState[i].U);
  }
  fvm_Release((void *)IO->InitState);
  fvm_Release((void **)IO->InitType);
  fvm_Release2((void **)IO->InitGroups);
  fvm_Release((void *)IO->nEGroups);
  
  fvm_Release((void *)IO);
  
  return 0;
}

/****************************************************************************/
//Function DestroyAll
int fvm_DestroyAll(fvm_All *All)
{
  fvm_DestroyData(All);
  fvm_DestroyMesh(All->Mesh);
  fvm_DestroyIO(All->IO);
  
  fvm_Release((void *)All);
  
  return 0;
}


