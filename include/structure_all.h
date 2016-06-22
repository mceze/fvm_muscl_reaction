/*
 *  structure_all.h
 *  fvm
 *
 *  Created by Marco Ceze on 7/22/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

/*This file declares the All structure*/

#ifndef _fvm_structure_all_h
#define _fvm_structure_all_h 1

#include "structure_mesh.h"
#include "structure_data.h"
#include "structure_IO.h"

/****************************************************************************/
//All Structure
typedef struct 
  {
    fvm_Mesh *Mesh;
    fvm_Data *Data;
    fvm_IO *IO;
  }
  fvm_All;


#endif