/*
 *  FluxER.c
 *  project1
 *
 *  Created by Marco Ceze on 2/14/09.
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

static int FluxER(double StateL[3], double StateR[3], double StateI[3]);

int fvm_FluxER(fvm_All *All, fvm_State StateL, fvm_State StateR, fvm_Face Face, double NormalFlux[4], double *MaxWaveSpeed)
{
  //This function computes the state at the interface using the 
  //Exact solution to the Riemann problem
  int i;
  double state_l[3], state_r[3], state_i[3], Normal[2], Tangent[2], rho, q, r, p;
  double F[4], G[4];
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  fvm_IO *IO;
  fvm_State InterfaceState;
  
  Data = All->Data;
  Mesh = All->Mesh;
  IO = All->IO;
  
  fvm_Alloc((void **)&InterfaceState.W, Mesh->Dim+2, sizeof(double));
  fvm_Alloc((void **)&InterfaceState.U, Mesh->Dim+2, sizeof(double));
  
  if (Face.Group != INTERIORFACE){
    for (i = 0; i < 2; i++){
      Normal[i] = Data->BFaceGroup[Face.Group].Normal[Face.Number][i];
      Tangent[i] = Data->BFaceGroup[Face.Group].Tangent[Face.Number][i];
    }
  }
  else {
    for (i = 0; i < 2; i++){
      Normal[i] = Data->IFace.Normal[Face.Number][i];
      Tangent[i] = Data->IFace.Tangent[Face.Number][i];
    }
  }
  
  state_l[0] = StateL.W[0];//rho
  state_l[1] = StateL.W[1]*Normal[0] + StateR.W[2]*Normal[1];//q
  state_l[2] = StateL.W[3];//p
  
  state_r[0] = StateR.W[0];//rho
  state_r[1] = StateR.W[1]*Normal[0] + StateR.W[2]*Normal[1];//q
  state_r[2] = StateR.W[3];//p
  
  FluxER(state_l, state_r, state_i);
  
  rho = state_i[0];
  q = state_i[1];
  p = state_i[2];
  if (q >= 0.0) r = StateL.W[1]*Tangent[0] + StateL.W[2]*Tangent[1];
  else r = StateR.W[1]*Tangent[0] + StateR.W[2]*Tangent[1];
  MaxWaveSpeed[0] = q + sqrt(GAMMA*p/rho);
  
  InterfaceState.W[0] = rho;
  InterfaceState.W[1] = q*Normal[0]+r*Tangent[0];
  InterfaceState.W[2] = q*Normal[1]+r*Tangent[1];
  InterfaceState.W[3] = p;
  fvm_ConvertState(InterfaceState, CONSERVED, Mesh->Dim, IO->nChemSpecies);
  
  fvm_CompFlux(InterfaceState, F, G);
  
  for (i = 0; i < 4; i++)
    NormalFlux[i] = F[i]*Normal[0] + G[i]*Normal[1];
  
  fvm_Release((void *)InterfaceState.W);
  fvm_Release((void *)InterfaceState.U);
  
  return 0;
}

static int FluxER(double StateL[3], double StateR[3], double StateI[3])
{
  //This function computes the state at the interface using the 
  //Exact solution to the Riemann problem
  int i;
  double uL, uR, pL, pR, rhoL, rhoR, aL, aR, mL, mR, p, u, rhoLs, rhoRs;
  double aLs, aRs, csL, csR, x1, x2, x3, x4, x5;
  double pUp, pLo, p_1, f_1, p1, rhs;
  
  //Unwrapping the states
  rhoL = StateL[0];
  uL = StateL[1];
  pL = StateL[2];
  
  rhoR = StateR[0];
  uR = StateR[1];
  pR = StateR[2];
  
  //computing the sound speeds
  aL = sqrt(GAMMA*pL/rhoL);
  aR = sqrt(GAMMA*pR/rhoR);
  
  //checking for equal states
  if(uL == uR && pL == pR && rhoL == rhoR){
    StateI[0] = rhoL;
    StateI[1] = uL;
    StateI[2] = pL;
    return 0;
  }
  
  //computing the sound speeds
  aL = sqrt(GAMMA*pL/rhoL);
  aR = sqrt(GAMMA*pR/rhoR);
  
  //Checking for cavitation
  if (uL + 2.0*aL/(GAMMA-1.0) < uR - 2.0*aR/(GAMMA-1.0)){
    printf("Cavitation detected! Exiting!\n");
    exit(1);
  }
  
  //Determining the pressure bounds
  //Upper limit given by a Poisson curve
  pUp = powl(((0.5*(GAMMA-1.0)*(uL-uR)+aL+aR)/(aL/(powl(pL,(0.5*(GAMMA-1.0)/GAMMA)))\
                                               + aR/(powl(pR,(0.5*(GAMMA-1.0)/GAMMA))))),(2.0*GAMMA/(GAMMA-1.0)));
  //Lower limit given by linear theory OR epsilon
  pLo = (aR*pL*rhoR+aL*rhoL*(pR+aR*rhoR*(uL-uR)))/(aL*rhoL+aR*rhoR);
  if(pLo < EPS) pLo = EPS;
  
  //initial estimate of pressure
  p = (pUp+pLo)/2.0; 
  
  p_1 = pUp;
  
  mL = fvm_MFlux(p_1, rhoL, aL, pL);
  mR = fvm_MFlux(p_1, rhoR, aR, pR);
  
  f_1 = fvm_SecRHS(p_1, mR, mL, uR, uL, pR, pL);
  
  
  for(i = 0; i < RSitmax; i++){
    rhs = fvm_SecRHS(p, mR, mL, uR, uL, pR, pL);
    p1 = p - rhs*(p - p_1)/(rhs - f_1);//secant method
    
    p_1 = p;
    f_1 = fvm_SecRHS(p, mR, mL, uR, uL, pR, pL);
    
    p = p1;
    
    mL = fvm_MFlux(p, rhoL, aL, pL);
    mR = fvm_MFlux(p, rhoR, aR, pR);
    
    rhs = fvm_SecRHS(p1, mR, mL, uR, uL, pR, pL);
    if (fabs(rhs) < EPS) break;
  }
  p = p1;
  //checking the final pressure
  if(p < 0.0){
    printf("Cavitation detected! Exiting!\n");
    exit(1);
  }
  
  u = (mL*uL + mR*uR - (pR - pL))/(mL + mR);
  
  rhoLs = fvm_Dens(p, rhoL, pL);
  rhoRs = fvm_Dens(p, rhoR, pR);
  aLs = sqrt(GAMMA*p/rhoLs);
  aRs = sqrt(GAMMA*p/rhoRs);
  
  if((p >= pL) && (p >= pR)){ //double compression
    csL = uL - aL*sqrt(1.0+((p-pL)/pL)*(GAMMA+1.0)/(2.0*GAMMA));
    csR = uR + aR*sqrt(1.0+((p-pR)/pR)*(GAMMA+1.0)/(2.0*GAMMA));
    x1 = csL;
    x2 = u;
    x3 = csR;
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0){
      StateI[0] = rhoR;
      StateI[1] = uR;
      StateI[2] = pR;
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 >= 0.0){
      StateI[0] = rhoRs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 <= 0.0 && x2 >= 0.0 && x3 >= 0.0){
      StateI[0] = rhoLs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 >= 0.0 && x2 >= 0.0 && x3 >= 0.0){
      StateI[0] = rhoL;
      StateI[1] = uL;
      StateI[2] = pL;
    }
    //printf("2comp\n");
  }
  if((p <= pL) && (p <= pR)){ //double expansion
    x1 = (uL - aL);
    x2 = (u - aLs);
    x3 = u;
    x4 = (u + aRs);
    x5 = (uR + aR);
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 <= 0.0 && x5 <= 0.0){
      StateI[0] = rhoR;
      StateI[1] = uR;
      StateI[2] = pR;
    }  
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 <= 0.0 && x5 >= 0.0){
      fvm_ExpFan(uR, aR, pR, rhoR, u, aRs, StateI, 1);
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 >= 0.0 && x5 >= 0.0){
      StateI[0] = rhoRs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 >= 0.0 && x4 >= 0.0 && x5 >= 0.0){
      StateI[0] = rhoLs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 <= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0 && x5 >= 0.0){
      fvm_ExpFan(uL, aL, pL, rhoL, u, aLs, StateI, -1);
    }
    if(x1 >= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0 && x5 >= 0.0){
      StateI[0] = rhoL;
      StateI[1] = uL;
      StateI[2] = pL;
    }
    //printf("2exp\n");
  }
  if((p <= pL) && (p >= pR)){ //expansion-compression
    csR = uR + aR*sqrt(1.0+((p-pR)/pR)*(GAMMA+1.0)/(2.0*GAMMA));
    x1 = (uL - aL);
    x2 = (u - aLs);
    x3 = u;
    x4 = csR;
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 <= 0.0){
      StateI[0] = rhoR;
      StateI[1] = uR;
      StateI[2] = pR;
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 >= 0.0){
      StateI[0] = rhoRs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      StateI[0] = rhoLs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 <= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      fvm_ExpFan(uL, aL, pL, rhoL, u, aLs, StateI, -1);
    }
    if(x1 >= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      StateI[0] = rhoL;
      StateI[1] = uL;
      StateI[2] = pL;
    }
    //printf("expcomp\n");
  }
  if((p >= pL) && (p <= pR)){ //compression-expansion
    csL = uL - aL*sqrt(1.0+((p-pL)/pL)*(GAMMA+1.0)/(2.0*GAMMA));
    x1 = csL;
    x2 = u;
    x3 = (u + aRs);
    x4 = (uR + aR);
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 <= 0.0){
      StateI[0] = rhoR;
      StateI[1] = uR;
      StateI[2] = pR;
    }    
    if(x1 <= 0.0 && x2 <= 0.0 && x3 <= 0.0 && x4 >= 0.0){
      fvm_ExpFan(uR, aR, pR, rhoR, u, aRs,StateI, 1);
    }
    if(x1 <= 0.0 && x2 <= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      StateI[0] = rhoRs;
      StateI[1] = u;
      StateI[2] = p;
    }    
    if(x1 <= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      StateI[0] = rhoLs;
      StateI[1] = u;
      StateI[2] = p;
    }
    if(x1 >= 0.0 && x2 >= 0.0 && x3 >= 0.0 && x4 >= 0.0){
      StateI[0] = rhoL;
      StateI[1] = uL;
      StateI[2] = pL;
    }
    //printf("compexp\n");
  }
  
  return 0;
}

