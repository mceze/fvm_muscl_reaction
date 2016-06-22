/*
 *  structure_IO.h
 *  fvm
 *
 *  Created by Marco Ceze and Chih-Kuang Kuan on 12/2/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

/*This file declares the IO data structure*/

#ifndef _fvm_structure_IO_h
#define _fvm_structure_IO_h 1

#include "main.h"

/****************************************************************************/
typedef struct 
  {
    int nInitStates;
    fvm_State *InitState;
    int **InitGroups;
    int *nEGroups;
    int *InitType;  //"Primitive" or "Conserved"
    char MeshFile[LONGLEN];
    double CFL;
    int itmax;
    double CalcTime;
    int ConstDt;
    int Limiter;
    int FluxType;
    int TimeMarching;
    int WriteEvery;
    double VRef;
    double LRef;
    double Kappa;
    //chemistry parameters
    int nChemSpecies;
    double *MolWeight;
    double *HeatOfFormation;
    double *StoichCoeff;
    double A_f;
    double b_f;
    double Ta_f;
    double A_b;
    double b_b;
    double Ta_b;
  }
  fvm_IO;

#endif