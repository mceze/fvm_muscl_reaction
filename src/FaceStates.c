/*
 *  FaceStates.c
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
//Function fvm_CezePhi
double fvm_CezePhi(double r)
{
  double phi;
  
  if (r >= 0.0){
    if (r < 1.0){
      phi = r;
    }
    else {
      phi = 4.0*pow(r,2.0)/(pow((r+1.0),2.0)+pow(r-1.0,2.0));
    }
  }
  else {
    phi = 0.0;
  }
  
  return phi;
}

/****************************************************************************/
//Function fvm_VanLeerPhi
double fvm_VanLeerPhi(double r)
{
  double phi;
  
  if (r >= 0.0)
    phi = (r+fabs(r))/(1.0+r);
  else
    phi = 0.0;
  
  return phi;
}

/****************************************************************************/
//Function fvm_KorenPhi
double fvm_KorenPhi(double r)
{
  double phi;
  
  if (r >= 0.0)
    phi = (2.0*pow(r,2.0)+r)/(2.0*pow(r,2.0)-r+2.0);
  else
    phi = 0.0;
  
  return phi;
}

/****************************************************************************/
//Function fvm_VanAlbadaPhi
double fvm_VanAlbadaPhi(double r)
{
  double phi;
  
  phi = (pow(r,2.0)+r)/(1.0+pow(r,2.0));
  
  return phi;
}

/****************************************************************************/
//Function fvm_SuperbeePhi
double fvm_SuperbeePhi(double r)
{
  double phi;
  
  phi = fmax(0.0,fmax(fmin(2.0*r,1.0),fmin(r,2.0)));
  
  return phi;
}

/****************************************************************************/
//Function fvm_MinModPhi
double fvm_MinModPhi(double r)
{
  double phi;
  
  if (r > 0.0){
    if (fabs(r) > 1.0) phi = 1.0;
    else phi = r;
  }
  else phi = 0.0;
  
  return phi;
}


/****************************************************************************/
//Function fvm_CalculatePhi
double fvm_CalculatePhi(fvm_IO *IO, double r)
{
  double phi;
  
  //printf("r: %lf\n",r);
  
  if (IO->Limiter == KOREN) phi = fvm_KorenPhi(r);
  if (IO->Limiter == VANALBADA) phi = fvm_VanAlbadaPhi(r);
  if (IO->Limiter == SUPERBEE) phi = fvm_SuperbeePhi(r);
  if (IO->Limiter == MINMOD) phi = fvm_MinModPhi(r);
  if (IO->Limiter == VANLEER) phi = fvm_VanLeerPhi(r);
  if (IO->Limiter == CEZE) phi = fvm_CezePhi(r);
  
  return phi;
}


/****************************************************************************/
//Function fvm_CompFaceStates
int fvm_CompFaceStates(fvm_All *All, fvm_Face Face, fvm_State StateL, fvm_State StateR)
{
  fvm_IO *IO;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  
  int egrpL, egrpR, elemL, elemR, faceL, faceR, i, k;
  int egrpLL, egrpRR, elemLL, elemRR, flag, oppface;
  int faceLL, faceRR;
  double rL, rR, phi;
  fvm_Face TempFace;
  fvm_State CellStateL, CellStateR, CellStateLL, CellStateRR;
  
  IO = All->IO;
  Mesh = All->Mesh;
  Data = All->Data;
  
  fvm_Alloc((void **)&CellStateL.W, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateL.U, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateLL.W, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateLL.U, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateR.W, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateR.U, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateRR.W, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&CellStateRR.U, Mesh->Dim+2+IO->nChemSpecies, sizeof(double));
  
  
  /****************************************************************************/
  //first order
  if (IO->Limiter == 0){
    if (Face.Group == INTERIORFACE){
      egrpL = Mesh->IFace[Face.Number].ElemGroupL;
      egrpR = Mesh->IFace[Face.Number].ElemGroupR;
      elemL = Mesh->IFace[Face.Number].ElemL;
      elemR = Mesh->IFace[Face.Number].ElemR;
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        StateL.W[k] = Data->ElemGroup[egrpL].State[elemL].W[k];
        StateR.W[k] = Data->ElemGroup[egrpR].State[elemR].W[k];
      }
    }
    else {
      egrpL = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].ElemGroup;
      elemL = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].Elem;
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        StateL.W[k] = Data->ElemGroup[egrpL].State[elemL].W[k];
        StateR.W[k] = Data->BFaceGroup[Face.Group].State[Face.Number].W[k];
      }
    }
  }
  /****************************************************************************/
  //With limiter
  if (IO->Limiter != 0){
    if (Face.Group != INTERIORFACE){//Boundary Face
      egrpL = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].ElemGroup;
      elemL = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].Elem;
      faceL = Mesh->BFaceGroup[Face.Group].BFace[Face.Number].Face;
      
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        CellStateL.W[k] = Data->ElemGroup[egrpL].State[elemL].W[k];
        CellStateR.W[k] = Data->BFaceGroup[Face.Group].State[Face.Number].W[k];
      }
      
      //getting CellStateLL
      oppface = fvm_FindOppositeFace(faceL, egrpL, All);
      fvm_NeighborAcrossFace(Mesh, egrpL, elemL, oppface, &egrpLL, &elemLL, &faceLL, &flag);
      if (flag == 1){//opposite face is a boundary
        TempFace.Group = Mesh->ElemGroup[egrpL].Face[elemL][oppface].Group;
        TempFace.Number = Mesh->ElemGroup[egrpL].Face[elemL][oppface].Number;
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateLL.W[k] = Data->BFaceGroup[TempFace.Group].State[TempFace.Number].W[k];
      }
      else
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateLL.W[k] = Data->ElemGroup[egrpLL].State[elemLL].W[k];      
      
      //Computing the Left State
      for (i = 0; i < Mesh->Dim+2+IO->nChemSpecies; i++){
        rL = (CellStateR.W[i]-CellStateL.W[i])/(CellStateL.W[i]-CellStateLL.W[i]+EPS)+EPS;
        phi = fvm_CalculatePhi(IO, rL);
        if (IO->Kappa < 1.0/3.0+EPS || IO->Kappa > 1.0/3.0-EPS)
          StateL.W[i] = CellStateL.W[i]+0.5*phi*(CellStateL.W[i]-CellStateLL.W[i]);
        else
          StateL.W[i] = CellStateL.W[i]+(phi/4.0)*((1.0-IO->Kappa)*(CellStateL.W[i]-CellStateLL.W[i])
                                                +(1.0+IO->Kappa)*(CellStateR.W[i]-CellStateL.W[i]));
      }
      //applying bcs to face state
      fvm_ApplyBC(Face, All, StateL, StateR);
    }
    else {//Interior Face
      elemL = Mesh->IFace[Face.Number].ElemL;
      egrpL = Mesh->IFace[Face.Number].ElemGroupL;
      faceL = Mesh->IFace[Face.Number].FaceL;
      
      elemR = Mesh->IFace[Face.Number].ElemR;
      egrpR = Mesh->IFace[Face.Number].ElemGroupR;
      faceR = Mesh->IFace[Face.Number].FaceR;
      
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        CellStateL.W[k] = Data->ElemGroup[egrpL].State[elemL].W[k];
        CellStateR.W[k] = Data->ElemGroup[egrpR].State[elemR].W[k];
      }
        
      //getting CellState LL
      oppface = fvm_FindOppositeFace(faceL, egrpL, All);
      fvm_NeighborAcrossFace(Mesh, egrpL, elemL, oppface, &egrpLL, &elemLL, &faceLL, &flag);
      if (flag == 1){//opposite face is a boundary
        TempFace.Group = Mesh->ElemGroup[egrpL].Face[elemL][oppface].Group;
        TempFace.Number = Mesh->ElemGroup[egrpL].Face[elemL][oppface].Number;
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateLL.W[k] = Data->BFaceGroup[TempFace.Group].State[TempFace.Number].W[k];
      }
      else
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateLL.W[k] = Data->ElemGroup[egrpLL].State[elemLL].W[k];
      
      //getting CellState RR
      oppface = fvm_FindOppositeFace(faceR, egrpR, All);
      fvm_NeighborAcrossFace(Mesh, egrpR, elemR, oppface, &egrpRR, &elemRR, &faceRR, &flag);
      if (flag == 1){//opposite face is a boundary
        TempFace.Group = Mesh->ElemGroup[egrpR].Face[elemR][oppface].Group;
        TempFace.Number = Mesh->ElemGroup[egrpR].Face[elemR][oppface].Number;
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateRR.W[k] = Data->BFaceGroup[TempFace.Group].State[TempFace.Number].W[k];
      }
      else
        for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
          CellStateRR.W[k] = Data->ElemGroup[egrpRR].State[elemRR].W[k];
      
      //Computing the Left and Right States
      for (i = 0; i < Mesh->Dim+2+IO->nChemSpecies; i++){
        rL = (CellStateR.W[i]-CellStateL.W[i])/(CellStateL.W[i]-CellStateLL.W[i]+EPS)+EPS;
        phi = fvm_CalculatePhi(IO, rL);
        if (IO->Kappa < 1.0/3.0+EPS || IO->Kappa > 1.0/3.0-EPS)
          StateL.W[i] = CellStateL.W[i]+0.5*phi*(CellStateL.W[i]-CellStateLL.W[i]);
        else
          StateL.W[i] = CellStateL.W[i]+(phi/4.0)*((1.0-IO->Kappa)*(CellStateL.W[i]-CellStateLL.W[i])
                                                   +(1.0+IO->Kappa)*(CellStateR.W[i]-CellStateL.W[i]));
        
        rR = (CellStateRR.W[i]-CellStateR.W[i])/(CellStateR.W[i]-CellStateL.W[i]+EPS)+EPS;
        phi = fvm_CalculatePhi(IO, 1.0/rR);
        if (IO->Kappa < 1.0/3.0+EPS || IO->Kappa > 1.0/3.0-EPS)
          StateR.W[i] = CellStateR.W[i]-0.5*rR*phi*(CellStateR.W[i]-CellStateL.W[i]);
        else
          StateR.W[i] = CellStateR.W[i]-(phi/4.0)*((1.0+IO->Kappa)*(CellStateR.W[i]-CellStateL.W[i])
                                                +(1.0-IO->Kappa)*(CellStateRR.W[i]-CellStateR.W[i]));
        
      }
    }
  }
  fvm_Release((void *)CellStateL.W);
  fvm_Release((void *)CellStateLL.W);
  fvm_Release((void *)CellStateR.W);
  fvm_Release((void *)CellStateRR.W);
  
  fvm_Release((void *)CellStateL.U);
  fvm_Release((void *)CellStateLL.U);
  fvm_Release((void *)CellStateR.U);
  fvm_Release((void *)CellStateRR.U);
  
  //Keeping the states coherent
  fvm_ConvertState(StateL, PRIMITIVE, Mesh->Dim, IO->nChemSpecies);
  fvm_ConvertState(StateR, PRIMITIVE, Mesh->Dim, IO->nChemSpecies);
  
  return 0;
}

