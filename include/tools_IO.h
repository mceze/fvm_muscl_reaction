/*
 *  tools_IO.h
 *  fvm
 *
 *  Created by Marco Ceze and Chih-Kuang Kuan on 12/2/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#ifndef _fvm_tools_IO_h
#define _fvm_tools_IO_h 1

#include "structure_mesh.h"
#include "structure_data.h"
#include "structure_IO.h"
#include "structure_all.h"

/****************************************************************************/
//Function fvm_ReadInput
int fvm_ReadInput(fvm_All  *All);

/****************************************************************************/
//Function fvm_ReadMesh
int fvm_ReadMesh(fvm_Mesh *Mesh, char MeshFile[]);

/****************************************************************************/
//Function fvm_WriteSolution
int fvm_WriteSolution(fvm_All *All, int it, double time);

#endif