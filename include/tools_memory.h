/*
 *  tools_memory.h
 *  fvm
 *
 *  Created by Marco Ceze on 7/22/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

/*This file declares the functions dealing with memory*/

#ifndef _fvm_tools_memory_h
#define _fvm_tools_memory_h 1

#include "structure_mesh.h"
#include "structure_data.h"
#include "structure_IO.h"
#include "structure_all.h"

/****************************************************************************/
int fvm_Alloc(void **pchunk, int n, int size);

/****************************************************************************/
int fvm_Alloc2( void ***pchunk, int n1, int n2, int size);

/****************************************************************************/
int fvm_CreateMesh(fvm_Mesh **Mesh);

/****************************************************************************/
int fvm_CreateData(fvm_Data **Data);

/****************************************************************************/
int fvm_CreateIO(fvm_IO **IO);

/****************************************************************************/
int fvm_CreateAll(fvm_All **pAll);

/****************************************************************************/
void fvm_Release(void *chunk);

/****************************************************************************/
void fvm_Release2(void **chunk);

/****************************************************************************/
void fvm_Release3(void ***chunk);

/****************************************************************************/
int fvm_DestroyMesh(fvm_Mesh *Mesh);

/****************************************************************************/
int fvm_DestroyData(fvm_All *All);

/****************************************************************************/
int fvm_DestroyIO(fvm_IO *IO);

/****************************************************************************/
int fvm_DestroyAll(fvm_All *All);

#endif

