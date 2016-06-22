/*
 *  structure_mesh.h
 *  fvm
 *
 *  Created by Marco Ceze on 7/20/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

/*This file declares the mesh structure*/

#ifndef _fvm_structure_mesh_h
#define _fvm_structure_mesh_h 1

#include "main.h"

/****************************************************************************/
//Face Structure
typedef struct
  {
    int Group;
    int Number;
  } 
  fvm_Face;

/****************************************************************************/
//Interior Face Structure
typedef struct
  {
    int ElemGroupL;
    int ElemL;
    int FaceL;
    int ElemGroupR;
    int ElemR;
    int FaceR;
    
    int OrientL;
    int OrientR;
  } 
  fvm_IFace;

/****************************************************************************/
//Boundary Face Structure
typedef struct
  {
    int ElemGroup;
    int Elem;
    int Face;
    
    int Orient;
  } 
  fvm_BFace;

/****************************************************************************/
//Boundary Face-group Structure
typedef struct
  {
    char Title[MEDIUMLEN];
    int nFace;
    int Type;
    fvm_BFace *BFace;
  } 
  fvm_BFaceGroup;

/****************************************************************************/
//Element-Group Structure
typedef struct 
  {
    int nElem;
    char Title[MEDIUMLEN];
    
    int nFace;
    fvm_Face **Face;
    
    int nNode;
    int **Node;
  } 
  fvm_ElemGroup;

/****************************************************************************/
//Mesh Structure
typedef struct
  {    
    int Dim;
    int nNode;
    
    double **Coord;
    
    int nIFace;
    fvm_IFace *IFace;
    
    int nBFaceGroup;
    fvm_BFaceGroup *BFaceGroup;
    
    int nElemGroup;
    fvm_ElemGroup *ElemGroup;
  }
  fvm_Mesh;

#endif
