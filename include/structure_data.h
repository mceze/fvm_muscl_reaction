/*
 *  structure_data.h
 *  fvm
 *
 *  Created by Marco Ceze on 7/21/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

/*This file declares the data structure*/

#ifndef _fvm_structure_data_h
#define _fvm_structure_data_h 1

/****************************************************************************/
//State-data Structure
typedef struct
  {
    double *U; //Array of conserved variables
    double *W; //Array of primitive variables
  } 
  fvm_State;

/****************************************************************************/
//Element-group-data Structure
typedef struct
  {
    double **ElemCenter; //ElemCenter[elem][dim]
    double *Dt;
    fvm_State *State;    //State[elem]
    double **dU_dt;
    double *Volume;
    //Jameson's-method-related data
    double **C;
    double **D;
    double *A;
    double *nu;
    double **Nabla2_Q;
  }
  fvm_EGrpData;

/****************************************************************************/
//Face-group-data Structure
typedef struct 
  {
    double *Area;     //Area[face] face is in the global numbering
    double **Normal;  //Normal[face][dim]
    double **Tangent;  //Tangent[face][dim]
    fvm_State *State; //State[face]
  } 
  fvm_BFGrpData;

/****************************************************************************/
//Interior-face-data Structure
typedef struct 
  {
    double *Area;     //Area[face] face is in the global numbering
    double **Normal;  //Normal[face][dim]
    double **Tangent;  //Tangent[face][dim]
  } 
  fvm_IFaceData;

/****************************************************************************/
//Data Structure
typedef struct
  {
    fvm_EGrpData *ElemGroup;
    fvm_BFGrpData *BFaceGroup;
    fvm_IFaceData IFace;
    double CurrentTime;
  } 
  fvm_Data;

#endif