/*
 *  FluxRoe.c
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
//Function fvm_FluxRoe (Still 2D)
int fvm_FluxRoe(fvm_All *All, fvm_State StateL, fvm_State StateR, fvm_Face Face, double NormalFlux[7], double *MaxWaveSpeed)
{
  //This function computes the state at the interface using
  //Roe's approximate solution to the Riemann problem
  int i, k;
  double uL, uR, vL, vR, pL, pR, rhoL, rhoR, aL, aR, HL, HR, qL, qR;
  double rR, rL, rHat, qHat, FluxL[4], FluxR[4];
  double w, rhoHat, uHat, vHat, HHat, aHat, RHat[4][4], FL[4], FR[4], GL[4], GR[4], del;
  double lambdaHat[4], dlambda[4], lambdaHatAbs[4], lambdaR[4], lambdaL[4], DeltaV[4], Dp, Dq, Dr, Drho;
  double Normal[2], Tangent[2], s1, s2;
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  fvm_IO *IO;
  
  Data = All->Data;
  Mesh = All->Mesh;
  IO = All->IO;
  
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
  
  //printf("Normal: %f %f Tangent: %f %f\n", Normal[0], Normal[1], Tangent[0], Tangent[1]);
  
  //Unwrapping the states
  rhoL = StateL.W[0];
  uL = StateL.W[1];
  vL = StateL.W[2];
  pL = StateL.W[3];
  
  rhoR = StateR.W[0];
  uR = StateR.W[1];
  vR = StateR.W[2];
  pR = StateR.W[3];
  
  qL = uL*Normal[0] + vL*Normal[1];
  qR = uR*Normal[0] + vR*Normal[1];
  
  rL = uL*Tangent[0] + vL*Tangent[1];
  rR = uR*Tangent[0] + vR*Tangent[1];
  
  //computing the sound speeds
  aL = fvm_SSpeed((fvm_State)StateL);
  aR = fvm_SSpeed((fvm_State)StateR);
  //computing the enthalpies
  HL = StateL.U[3]/rhoL + pL/rhoL;
  HR = StateR.U[3]/rhoR + pR/rhoR;
  
  w = sqrt(rhoR/rhoL);
  //computing hat-states
  rhoHat = w*rhoL;
  uHat = (uL + w*uR)/(1.0 + w);
  vHat = (vL + w*vR)/(1.0 + w);
  HHat = (HL + w*HR)/(1.0 + w);
  qHat = uHat*Normal[0] + vHat*Normal[1];
  rHat = uHat*Tangent[0] + vHat*Tangent[1];
  aHat = sqrt((GAMMA-1.0)*(HHat-0.5*(pow(uHat,2.0)+pow(vHat,2.0))));
  MaxWaveSpeed[0] = fabs(qHat) + aHat;
  //right hat-eigenvectors matrix
  RHat[0][0] = 1.0;
  RHat[0][1] = 1.0;
  RHat[0][2] = 0.0;
  RHat[0][3] = 1.0;
  
  RHat[1][0] = uHat - aHat*Normal[0];
  RHat[1][1] = uHat;
  RHat[1][2] = -Normal[1];
  RHat[1][3] = uHat + aHat*Normal[0];
  
  RHat[2][0] = vHat + aHat*(-Normal[1]);
  RHat[2][1] = vHat;
  RHat[2][2] = Normal[0];
  RHat[2][3] = vHat - aHat*(-Normal[1]);
  
  RHat[3][0] = HHat - qHat*aHat;
  RHat[3][1] = 0.5*(pow(uHat,2.0)+pow(vHat,2.0));
  RHat[3][2] = rHat;
  RHat[3][3] = HHat + qHat*aHat;
  
  //lambdas
  lambdaL[0] = qL - aL;
  lambdaL[1] = qL;
  lambdaL[2] = qL;
  lambdaL[3] = qL + aL;
  
  lambdaR[0] = qR - aR;
  lambdaR[1] = qR;
  lambdaR[2] = qR;
  lambdaR[3] = qR + aR;
  
  lambdaHat[0] = qHat - aHat;
  lambdaHat[1] = qHat;
  lambdaHat[2] = qHat;
  lambdaHat[3] = qHat + aHat;
  
  dlambda[0] = 2.0*fmin(aHat,fmax(0.0, 2.0*(lambdaR[0] - lambdaL[0])));
  dlambda[1] = 0.0;
  dlambda[2] = 0.0;
  dlambda[3] = 2.0*fmin(aHat,fmax(0.0, 2.0*(lambdaR[3] - lambdaL[3])));
  
  /* for(i = 0; i < 4; i++){
      if(fabs(lambdaHat[i]) >= dlambda[i]/2.0) lambdaHatAbs[i] = fabs(lambdaHat[i]);
      else lambdaHatAbs[i] = pow(lambdaHat[i],2.0)/dlambda[i] + dlambda[i]/4.0;
    } */
  for (i = 0; i < 4; i++){
    del = fmax(0.0, fmax((lambdaHat[i]-lambdaL[i]), (lambdaR[i]-lambdaHat[i])));
    
    if (fabs(lambdaHat[i])>=del)
      lambdaHatAbs[i] = fabs(lambdaHat[i]);
    else
      lambdaHatAbs[i] = (pow(lambdaHat[i],2.0)+pow(del,2.0))/(2.0*del);
  }
  
  //Delta V
  Dp = pR-pL;
  Dq = qR-qL;
  Dr = rR-rL;
  Drho = rhoR-rhoL;
  
  DeltaV[0] = (Dp-rhoHat*aHat*Dq)/(2.0*pow(aHat,2.0));
  DeltaV[1] = -(Dp-pow(aHat,2.0)*Drho)/pow(aHat,2.0);
  DeltaV[2] = rhoHat*Dr;
  DeltaV[3] = (Dp+rhoHat*aHat*Dq)/(2.0*pow(aHat,2.0));
  
  //computing the left and right fluxes
  fvm_CompFlux((fvm_State)StateL, FL, GL);
  fvm_CompFlux((fvm_State)StateR, FR, GR);
  for (i = 0; i < 4; i++){
    FluxL[i] = FL[i]*Normal[0] + GL[i]*Normal[1];
    FluxR[i] = FR[i]*Normal[0] + GR[i]*Normal[1];
  }
  //printf("FluxL: %f %f %f %f\n",FluxL[0],FluxL[1],FluxL[2],FluxL[3]);
  //printf("FluxR: %f %f %f %f\n",FluxR[0],FluxR[1],FluxR[2],FluxR[3]);
  
  //computing the interface flux vector
  if (All->IO->FluxType == ROEE){
    for(i = 0; i < 4; i++){
      NormalFlux[i] = 0.5*(FluxL[i]+FluxR[i]);
      for(k = 0; k < 4; k++){
        NormalFlux[i] += -0.5*lambdaHatAbs[k]*DeltaV[k]*RHat[i][k];
      }
    }
    //Chemical Species convection
    if (qHat >= 0)
      for (i = 0; i < IO->nChemSpecies; i++)
        NormalFlux[4+i] = qHat*StateL.U[4+i];
    else
      for (i = 0; i < IO->nChemSpecies; i++)
        NormalFlux[4+i] = qHat*StateR.U[4+i];
  }
  //printf("normal flux: %f %f %f %f %f %f %f\n",NormalFlux[0],NormalFlux[1],NormalFlux[2],NormalFlux[3],NormalFlux[4],NormalFlux[5],NormalFlux[6]);
  if (All->IO->FluxType == ROE){
    if (qHat < 0) s1 = -1.0;
    else s1 = 1.0;
    if ((pow(qHat,2.0)-pow(aHat,2.0)) < 0) s2 = -1.0;
    else s2 = 1.0;
    NormalFlux[0] = 0.5*(1.0+s1)*FluxL[0]+0.5*(1.0-s1)*FluxR[0]-0.5*(1.0-s2)*(qHat-s1*aHat)*((rhoHat*aHat*Dq-s1*Dp)/(2.0*pow(aHat, 2.0)))*(1.0);
    NormalFlux[1] = 0.5*(1.0+s1)*FluxL[1]+0.5*(1.0-s1)*FluxR[1]-0.5*(1.0-s2)*(qHat-s1*aHat)*((rhoHat*aHat*Dq-s1*Dp)/(2.0*pow(aHat, 2.0)))*(qHat-s1*aHat);
    NormalFlux[2] = 0.5*(1.0+s1)*FluxL[2]+0.5*(1.0-s1)*FluxR[2]-0.5*(1.0-s2)*(qHat-s1*aHat)*((rhoHat*aHat*Dq-s1*Dp)/(2.0*pow(aHat, 2.0)))*(qHat-s1*aHat);
    NormalFlux[3] = 0.5*(1.0+s1)*FluxL[3]+0.5*(1.0-s1)*FluxR[3]-0.5*(1.0-s2)*(qHat-s1*aHat)*((rhoHat*aHat*Dq-s1*Dp)/(2.0*pow(aHat, 2.0)))*(HHat-s1*qHat*aHat);
    //Chemical Species convection
    if (qHat >= 0)
      for (i = 0; i < IO->nChemSpecies; i++)
        NormalFlux[4+i] = qHat*StateL.U[4+i];
    else
      for (i = 0; i < IO->nChemSpecies; i++)
        NormalFlux[4+i] = qHat*StateR.U[4+i];
  
  }
  //printf("Normal: %f %f\n", Normal[0], Normal[1]);
  //printf("%f %f %f %f\n",NormalFlux[0],NormalFlux[1],NormalFlux[2],NormalFlux[3]);
  return 0;
}

