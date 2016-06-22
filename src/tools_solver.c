/*
 *  tools_solver.c
 *  fvm
 *
 *  Created by Marco Ceze on 8/3/09.
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
//Function fvm_ConvertState
int fvm_ConvertState(fvm_State State, int Source, int Dim, int nChemSpecies)
{
  int i;
  double sum;
  
  if (Dim == 2){
    if (Source == CONSERVED){
      State.W[0] = State.U[0];
      State.W[1] = State.U[1]/State.U[0];
      State.W[2] = State.U[2]/State.U[0];
      State.W[3] = (GAMMA-1.0)*(State.U[3]-State.W[0]*
                                (pow(State.W[1],2.0)+pow(State.W[2],2.0))/2.0);
      for (i = 0; i < nChemSpecies; i++)
        State.W[4+i] = State.U[4+i]/State.U[0];
    }
    else if (Source == PRIMITIVE){
      //if (State.W[0] < EPS) State.W[0] = EPS;
      //if (State.W[3] < EPS) State.W[3] = EPS;
      State.U[0] = State.W[0];
      State.U[1] = State.U[0]*State.W[1];
      State.U[2] = State.U[0]*State.W[2];
      State.U[3] = State.W[3]/(GAMMA-1.0)+State.W[0]*(pow(State.W[1],2.0)
                                                      +pow(State.W[2],2.0))/2.0;
      for (i = 0; i < nChemSpecies; i++)
        State.U[4+i] = State.W[4+i]*State.W[0];
    }
  }
  else{
    printf("Dimension not supported! Exiting!");
    exit(0);
  }
  //normalizing mass fractions
  sum = 0.0;
  for (i = 0; i < nChemSpecies; i++)
    sum += State.W[4+i];
  for (i = 0; i < nChemSpecies; i++)
    State.W[4+i] = State.W[4+i]/sum; 
  
  //normalizing species densities
  sum = 0.0;
  for (i = 0; i < nChemSpecies; i++)
    sum += State.U[4+i];
  for (i = 0; i < nChemSpecies; i++)
    State.U[4+i] = State.U[0]*State.U[4+i]/sum;
  
  return 0;
}

/****************************************************************************/
//Function fvm_SSpeed
double fvm_SSpeed(fvm_State State)
{
  //This function computes the sound speed for p and rho (2D!!!!)
  double a, p, rho;
  
  rho = State.W[0];
  p = State.W[3];
  a = sqrt(GAMMA*p/rho);
  
  return a;
}

/****************************************************************************/
//Function fvm_CompFlux
int fvm_CompFlux(fvm_State State, double F[4], double G[4])
{
  double rho, u, v, p, H, E;
  
  rho = State.W[0];
  u = State.W[1];
  v = State.W[2];
  p = State.W[3];
  E = State.U[3]/rho;
  H = E + p/rho;
  
  F[0] = rho*u;
  F[1] = rho*pow(u,2.0) + p;
  F[2] = rho*u*v;
  F[3] = rho*u*H;
  
  G[0] = rho*v;
  G[1] = rho*u*v;
  G[2] = rho*pow(v,2.0) + p;
  G[3] = rho*v*H;
  
  return 0;
}

/****************************************************************************/
//Function fvm_ApplyBC
int fvm_ApplyBC(fvm_Face Face, fvm_All *All, fvm_State InsideState, fvm_State OutsideState)
{
  int i, fgrp, face;
  double Normal[All->Mesh->Dim], Tangent[All->Mesh->Dim], qin, rin, qout, rout;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  fvm_IO *IO;
  
  Mesh = All->Mesh;
  Data = All->Data;
  IO = All->IO;
  
  fgrp = Face.Group;
  face = Face.Number;
  for (i = 0; i < Mesh->Dim; i++){
    Normal[i] = Data->BFaceGroup[fgrp].Normal[face][i];
    Tangent[i] = Data->BFaceGroup[fgrp].Tangent[face][i];
  }
  //This is specific for 2D
  qin = InsideState.W[1]*Normal[0]+InsideState.W[2]*Normal[1];
  rin = InsideState.W[1]*Tangent[0]+InsideState.W[2]*Tangent[1];
  
  if (Mesh->BFaceGroup[fgrp].Type == WALL){
    qout = -qin;
    rout = rin;
    OutsideState.W[0] = InsideState.W[0];
    OutsideState.W[1] = qout*Normal[0] + rout*Tangent[0];
    OutsideState.W[2] = qout*Normal[1] + rout*Tangent[1];    
    OutsideState.W[3] = InsideState.W[3];
    for (i = 0; i < IO->nChemSpecies; i++)
      OutsideState.W[4+i] = InsideState.W[4+i];
  }
  if (Mesh->BFaceGroup[fgrp].Type == INFLOWOUTFLOW){
    for (i = 0; i < (Mesh->Dim+2)+IO->nChemSpecies; i++){
      OutsideState.W[i] = InsideState.W[i];
    }
  }
  if (Mesh->BFaceGroup[fgrp].Type == FARFIELD){
    printf("BC not implemented yet! Exiting!\n");
    exit(0);
  }
  
  fvm_ConvertState(OutsideState, PRIMITIVE, Mesh->Dim, IO->nChemSpecies);
  
  return 0;
}

/****************************************************************************/
//Function fvm_BCond
int fvm_BCond(fvm_All *All)
{
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  int fgrp, face, elemL, egrpL, faceL;
  
  Mesh = All->Mesh;
  Data = All->Data;
  
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
      egrpL = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
      elemL = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
      faceL = Mesh->BFaceGroup[fgrp].BFace[face].Face;
      fvm_ApplyBC(Mesh->ElemGroup[egrpL].Face[elemL][faceL], All, 
                  Data->ElemGroup[egrpL].State[elemL], 
                  Data->BFaceGroup[fgrp].State[face]);
    }
  }
  
  return 0;
}

/****************************************************************************/
//Function fvm_CompEulerRHS
int fvm_CompEulerRHS(fvm_All *All)
{
  int egrpL, elemL, faceL, egrpR, elemR, faceR, fgrp;
  int egrp, elem, face, k;
  double *NormalFlux, MaxWaveSpeed, **AveMaxWaveSpeed;
  fvm_Face Face;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  fvm_IO *IO;
  fvm_State StateL, StateR;
  
  Mesh = All->Mesh;
  Data = All->Data;
  IO = All->IO;
  
  fvm_Alloc((void **)&StateL.W, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
  fvm_Alloc((void **)&StateL.U, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
  fvm_Alloc((void **)&StateR.W, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
  fvm_Alloc((void **)&StateR.U, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
  fvm_Alloc((void **)&NormalFlux, Mesh->Dim+2+IO->nChemSpecies,sizeof(double));
  
  //updating boundary conditions in ghost cells
  //this is important for the face value calculation
  fvm_BCond(All);
  
  
  AveMaxWaveSpeed = malloc(Mesh->nElemGroup * sizeof(double));
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    AveMaxWaveSpeed[egrp] = malloc(Mesh->ElemGroup[egrp].nElem * sizeof(double));
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      AveMaxWaveSpeed[egrp][elem] = 0.0;
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        Data->ElemGroup[egrp].dU_dt[elem][k] = 0.0;
      }
    }
  }
  
  //Boundary faces contribution
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
      faceL = Mesh->BFaceGroup[fgrp].BFace[face].Face;
      egrpL = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
      elemL = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
      
      MaxWaveSpeed = 0.0;
      
      Face.Group = Mesh->ElemGroup[egrpL].Face[elemL][faceL].Group;
      Face.Number = Mesh->ElemGroup[egrpL].Face[elemL][faceL].Number;
      
      fvm_CompFaceStates(All, Face, StateL, StateR);
      if (IO->FluxType == ROE || IO->FluxType == ROEE) 
        fvm_FluxRoe(All, StateL, StateR, Face, NormalFlux, &MaxWaveSpeed);
      else if(IO->FluxType == ER)
        fvm_FluxER(All, StateL, StateR, Face, NormalFlux, &MaxWaveSpeed);
      //printf("NormalFlux = %f %f %f %f\n", NormalFlux[0], NormalFlux[1], NormalFlux[2], NormalFlux[3]);
      
      for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
        Data->ElemGroup[egrpL].dU_dt[elemL][k] += -NormalFlux[k]*Data->BFaceGroup[fgrp].Area[face]
        /Data->ElemGroup[egrpL].Volume[elemL];
        //printf("%f\n",Data->ElemGroup[egrpL].dU_dt[elemL][k]);
      }
      AveMaxWaveSpeed[egrpL][elemL] += 0.5*MaxWaveSpeed*Data->BFaceGroup[fgrp].Area[face]/Data->ElemGroup[egrpL].Volume[elemL];
    }
  }
  
  //Interior faces contribution
  for (face = 0; face < Mesh->nIFace; face++){
    faceL = Mesh->IFace[face].FaceL;
    egrpL = Mesh->IFace[face].ElemGroupL;
    elemL = Mesh->IFace[face].ElemL;
    
    faceR = Mesh->IFace[face].FaceR;
    egrpR = Mesh->IFace[face].ElemGroupR;
    elemR = Mesh->IFace[face].ElemR;
    
    Face.Group = Mesh->ElemGroup[egrpL].Face[elemL][faceL].Group;
    Face.Number = Mesh->ElemGroup[egrpL].Face[elemL][faceL].Number;
    
    fvm_CompFaceStates(All, Face, StateL, StateR);
    if (IO->FluxType == ROE || IO->FluxType == ROEE) 
      fvm_FluxRoe(All, StateL, StateR, Face, NormalFlux, &MaxWaveSpeed);
    else if(IO->FluxType == ER)
      fvm_FluxER(All, StateL, StateR, Face, NormalFlux, &MaxWaveSpeed);
    
    for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
      Data->ElemGroup[egrpL].dU_dt[elemL][k] += -NormalFlux[k]*Data->IFace.Area[face]
      /Data->ElemGroup[egrpL].Volume[elemL];
      
      Data->ElemGroup[egrpR].dU_dt[elemR][k] += NormalFlux[k]*Data->IFace.Area[face]
      /Data->ElemGroup[egrpR].Volume[elemR];
    }
    AveMaxWaveSpeed[egrpL][elemL] += 0.5*MaxWaveSpeed*Data->IFace.Area[face]/Data->ElemGroup[egrpL].Volume[elemL];
    AveMaxWaveSpeed[egrpR][elemR] += 0.5*MaxWaveSpeed*Data->IFace.Area[face]/Data->ElemGroup[egrpR].Volume[elemR];
    
    
  }
  
  //computing dt for each element
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++)
      Data->ElemGroup[egrp].Dt[elem] = All->IO->CFL/AveMaxWaveSpeed[egrp][elem];
    free(AveMaxWaveSpeed[egrp]);
  }
  free(AveMaxWaveSpeed);
  
  fvm_Release((void *)StateL.W);
  fvm_Release((void *)StateL.U);
  fvm_Release((void *)StateR.W);
  fvm_Release((void *)StateR.U);
  fvm_Release((void *)NormalFlux);
  
  return 0;
}

/****************************************************************************/
//Function fvm_CompChemRHS
int fvm_CompChemRHS(fvm_All *All)
{
  //This function assumes 3 species!!!
  int egrp, elem, i;
  double T, kf, kb, C[3], dY_dt[3], Dt_chem, dt, dQ_dt;
  fvm_Mesh *Mesh;
  fvm_IO *IO;
  fvm_Data *Data;
  
  fvm_Alloc((void **)&C, IO->nChemSpecies, sizeof(double));
  fvm_Alloc((void **)&dY_dt, IO->nChemSpecies, sizeof(double));
  
  Mesh = All->Mesh;
  IO = All->IO;
  Data = All->Data;
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      //Non-dimensional temperature (T = gamma*P/rho)
      T = GAMMA*Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+1]/Data->ElemGroup[egrp].State[elem].W[0];
      //rate coefficients (modified Arrhenius Law)
      kf = IO->A_f*pow(T,IO->b_f)*exp(-IO->Ta_f/T);
      kb = IO->A_b*pow(T,IO->b_b)*exp(-IO->Ta_b/T);
      for (i = 0; i < IO->nChemSpecies; i++)
        C[i] = Data->ElemGroup[egrp].State[elem].U[Mesh->Dim+2+i]/IO->MolWeight[i];
      
      //Assuming only 3 species
      // nu_A*A + nu_B*B <=> nu_C*C
      dY_dt[0] = IO->MolWeight[0]*IO->StoichCoeff[0]*(-kf*pow(C[0],IO->StoichCoeff[0])*pow(C[1],IO->StoichCoeff[1])+kb*pow(C[2],IO->StoichCoeff[2]));
      dY_dt[1] = IO->MolWeight[1]*IO->StoichCoeff[1]*(-kf*pow(C[0],IO->StoichCoeff[0])*pow(C[1],IO->StoichCoeff[1])+kb*pow(C[2],IO->StoichCoeff[2]));
      dY_dt[2] = -IO->MolWeight[2]*IO->StoichCoeff[2]*(-kf*pow(C[0],IO->StoichCoeff[0])*pow(C[1],IO->StoichCoeff[1])+kb*pow(C[2],IO->StoichCoeff[2]));
      
      //computing chemistry timestep
      Dt_chem = 1e6;
      for (i = 0; i < IO->nChemSpecies; i++){
        dt = Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+2+i]/fabs(dY_dt[i]);
        if (Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+2+i] < EPS)
          dt = EPS;
        if (dt < Dt_chem)
          Dt_chem = dt;
        
        Data->ElemGroup[egrp].dU_dt[elem][Mesh->Dim+2+i] += Data->ElemGroup[egrp].State[elem].U[0]*dY_dt[i];
        Data->ElemGroup[egrp].dU_dt[elem][0] += Data->ElemGroup[egrp].State[elem].U[0]*dY_dt[i];
      }
      //In the case where the chemistry dt is smaller we set the Dt at that cell as Dt_chem
      if (Data->ElemGroup[egrp].Dt[elem] > Dt_chem){
        Data->ElemGroup[egrp].Dt[elem] = Dt_chem;
      }
      //Heat of combustion
      dQ_dt = IO->HeatOfFormation[2]*Data->ElemGroup[egrp].State[elem].U[0]*dY_dt[2]/IO->MolWeight[2]-
              IO->HeatOfFormation[0]*Data->ElemGroup[egrp].State[elem].U[0]*dY_dt[0]/IO->MolWeight[0]-
              IO->HeatOfFormation[1]*Data->ElemGroup[egrp].State[elem].U[0]*dY_dt[1]/IO->MolWeight[1];
   
      Data->ElemGroup[egrp].dU_dt[elem][Mesh->Dim+1] += -dQ_dt;
      
    }
  }

  return 0;
}
/****************************************************************************/
//Function fvm_TimeMarching
int fvm_TimeMarching(fvm_All *All)
{
  int it, egrp, elem, k, write_flag;
  double ***Uo, Dt;
  fvm_IO *IO;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  
  IO = All->IO;
  Mesh = All->Mesh;
  Data = All->Data;
  
  //allocating Uo  
  Uo = malloc(Mesh->nElemGroup * sizeof(double));
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    Uo[egrp] = malloc(Mesh->ElemGroup[egrp].nElem * sizeof(double));
    for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++)
      Uo[egrp][elem] = malloc((Mesh->Dim+2+IO->nChemSpecies) * sizeof(double));
  }
  
  write_flag = 0;
  if (IO->FluxType == JAMESON){
    printf("FluxType not supported yet.\n");
    exit(0);
  }
  else if (IO->TimeMarching == RK3){
    it =0;
    fvm_WriteSolution(All, it, Data->CurrentTime);
    while (Data->CurrentTime < IO->CalcTime){
      
      //backing up state
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
          for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
            Uo[egrp][elem][k] = Data->ElemGroup[egrp].State[elem].U[k];
        }
      }
      //Computing RHS(U0)
      fvm_CompEulerRHS(All);
      fvm_CompChemRHS(All);
      fvm_CalcStepSize(All);
      Dt = Data->ElemGroup[0].Dt[0];
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
          for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
            Data->ElemGroup[egrp].State[elem].U[k] = Uo[egrp][elem][k] + (Dt/3.0)*Data->ElemGroup[egrp].dU_dt[elem][k];
          fvm_ConvertState(Data->ElemGroup[egrp].State[elem], CONSERVED, Mesh->Dim, IO->nChemSpecies);
        }
      }
      //Computing RHS(U1)
      fvm_CompEulerRHS(All);
      fvm_CompChemRHS(All);
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
          for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
            Data->ElemGroup[egrp].State[elem].U[k] = Uo[egrp][elem][k] + (Dt/2.0)*Data->ElemGroup[egrp].dU_dt[elem][k];
          fvm_ConvertState(Data->ElemGroup[egrp].State[elem], CONSERVED, Mesh->Dim, IO->nChemSpecies);
        }
      }
      //Computing RHS(U2)
      fvm_CompEulerRHS(All);
      fvm_CompChemRHS(All);
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
          for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
            Data->ElemGroup[egrp].State[elem].U[k] = Uo[egrp][elem][k] + Dt*Data->ElemGroup[egrp].dU_dt[elem][k];
          fvm_ConvertState(Data->ElemGroup[egrp].State[elem], CONSERVED, Mesh->Dim, IO->nChemSpecies);
        }
      }
            
      if (it/IO->WriteEvery >= write_flag){
        write_flag++;
        fvm_WriteSolution(All, it+1, Data->CurrentTime);
      }
      it++;
      printf("Dt=%e @ iteration: %d Current Time: %e\n", Dt, it, Data->CurrentTime);
    }
  }
  else if (IO->TimeMarching == FE){
    it =0;
    fvm_WriteSolution(All, it, Data->CurrentTime);
    while (Data->CurrentTime < IO->CalcTime){
      //printf("-------------------------------------------------------------\n");
      fvm_CompEulerRHS(All);
      //fvm_CompChemRHS(All);
      fvm_CalcStepSize(All);
      Dt = Data->ElemGroup[0].Dt[0];
      for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
        for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
          for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++)
            Data->ElemGroup[egrp].State[elem].U[k] += Dt*Data->ElemGroup[egrp].dU_dt[elem][k];
          fvm_ConvertState(Data->ElemGroup[egrp].State[elem], CONSERVED, Mesh->Dim, IO->nChemSpecies);
          //printf("Conserved State: %f %f %f %f\n",Data->ElemGroup[egrp].State[elem].U[0],Data->ElemGroup[egrp].State[elem].U[1],Data->ElemGroup[egrp].State[elem].U[2],Data->ElemGroup[egrp].State[elem].U[3]);
          //printf("Primitive State: %f %f %f %f\n",Data->ElemGroup[egrp].State[elem].W[0],Data->ElemGroup[egrp].State[elem].W[1],Data->ElemGroup[egrp].State[elem].W[2],Data->ElemGroup[egrp].State[elem].W[3]);
        }
      }
      if (it/IO->WriteEvery >= write_flag){
        write_flag++;
        fvm_WriteSolution(All, it+1, Data->CurrentTime);
      }
      it++;
      printf("Dt=%e @ iteration: %d Current Time: %e\n", Dt, it, Data->CurrentTime);
    }
  }
  
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++)
      free(Uo[egrp][elem]);
    free(Uo[egrp]);
  }
  free(Uo);
  
  fvm_WriteSolution(All, it, Data->CurrentTime);
  
  return 0;
}

/****************************************************************************/
//Function fvm_CalcStepSize
int fvm_CalcStepSize(fvm_All *All)
{
  int egrp, elem, face, f, fgrp;
  double MinDt, smallest_length, length, c, u, v;
  fvm_IO *IO;
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  
  IO = All->IO;
  Data = All->Data;
  Mesh = All->Mesh;
  
  //Calculating timestep size for JAMESON Flux
  if (IO->FluxType == JAMESON){
    for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
      for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
        smallest_length = 1.0e6;
        for (f = 0; f < Mesh->ElemGroup[egrp].nFace; f++){
          fgrp = Mesh->ElemGroup[egrp].Face[elem][f].Group;
          face = Mesh->ElemGroup[egrp].Face[elem][f].Number;
          if (fgrp == INTERIORFACE)
            length = Data->IFace.Area[face];
          else
            length = Data->BFaceGroup[fgrp].Area[face];
          if (length <= smallest_length)
            smallest_length = length;
        }
        u = Data->ElemGroup[egrp].State[elem].W[1];
        v = Data->ElemGroup[egrp].State[elem].W[2];
        c = sqrt(pow(u,2.0)+pow(v,2.0)) + fvm_SSpeed(Data->ElemGroup[egrp].State[elem]);
        Data->ElemGroup[egrp].Dt[elem] = smallest_length*(IO->CFL)/c;
      }
    }
  }
  
  //Homogenizing the time-step if requested
  if (IO->ConstDt == 1){
    //Temp change
    MinDt = 1.0e6;
    for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
      for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
        if (Data->ElemGroup[egrp].Dt[elem] < MinDt)
          MinDt = Data->ElemGroup[egrp].Dt[elem];
      }
    }        
    if(Data->CurrentTime + MinDt > IO->CalcTime){
      MinDt = IO->CalcTime - Data->CurrentTime;
    }
    Data->CurrentTime += MinDt;
    //printf("Dt: %f\n",MinDt);
    for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
      for(elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++)
        Data->ElemGroup[egrp].Dt[elem] = MinDt;
  }
  return 0;
}

/****************************************************************************/
//Function fvm_Convective
int fvm_Convective(fvm_All *All)
{
  int fgrp, face, egrpL, elemL, egrpR, elemR, k;
  double F[4], G[4];
  fvm_Data *Data;
  fvm_Mesh *Mesh;
  fvm_State InterfaceState;
  
  Mesh = All->Mesh;
  Data = All->Data;
  
  fvm_Alloc((void **)&InterfaceState.W, Mesh->Dim+2,sizeof(double));
  fvm_Alloc((void **)&InterfaceState.U, Mesh->Dim+2,sizeof(double));
  
  //BFaces
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
      egrpL = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
      elemL = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
      for (k = 0; k < Mesh->Dim+2; k++){
        InterfaceState.W[k] = 0.5*(Data->ElemGroup[egrpL].State[elemL].W[k]+
                                   Data->BFaceGroup[fgrp].State[face].W[k]);
        InterfaceState.U[k] = 0.5*(Data->ElemGroup[egrpL].State[elemL].U[k]+
                                   Data->BFaceGroup[fgrp].State[face].U[k]);
      }
      fvm_CompFlux(InterfaceState, F, G);
      
      for (k = 0; k < Mesh->Dim+2; k++)
        Data->ElemGroup[egrpL].C[elemL][k] = (F[k]*Data->BFaceGroup[fgrp].Normal[face][0]+
                                               G[k]*Data->BFaceGroup[fgrp].Normal[face][1])*
                                               Data->BFaceGroup[fgrp].Area[face];
    }
  }
  //IFaces
  for (face = 0; face < Mesh->nIFace; face++){
    egrpL = Mesh->IFace[face].ElemGroupL;
    elemL = Mesh->IFace[face].ElemL;
    egrpR = Mesh->IFace[face].ElemGroupR;
    elemR = Mesh->IFace[face].ElemR;
    
    for (k = 0; k < Mesh->Dim+2; k++){
      InterfaceState.W[k] = 0.5*(Data->ElemGroup[egrpL].State[elemL].W[k]+
                                 Data->ElemGroup[egrpR].State[elemR].W[k]);
      InterfaceState.U[k] = 0.5*(Data->ElemGroup[egrpL].State[elemL].U[k]+
                                 Data->ElemGroup[egrpR].State[elemR].U[k]);
    }
    fvm_CompFlux(InterfaceState, F, G);
    
    for (k = 0; k < Mesh->Dim+2; k++){
      Data->ElemGroup[egrpL].C[elemL][k] = (F[k]*Data->IFace.Normal[face][0]+
                                             G[k]*Data->IFace.Normal[face][1])*
                                             Data->BFaceGroup[fgrp].Area[face];
      Data->ElemGroup[egrpR].C[elemR][k] = -(F[k]*Data->IFace.Normal[face][0]+
                                             G[k]*Data->IFace.Normal[face][1])*
                                             Data->BFaceGroup[fgrp].Area[face];
    }
    
  }
  
  fvm_Release((void *)InterfaceState.W);
  fvm_Release((void *)InterfaceState.U);
   
  return 0;
}

/****************************************************************************/
//Function fvm_Dissipative
int fvm_Dissipative(fvm_All *All)
{
  int f, k, face, fgrp, n, neighbor, neighbor_grp, neighbor_face, egrp, elem, flag;
  double k2, k4, E2, E4, Normal[2], nu_neighbor, A_neighbor, Nabla2_Q_neighbor[4], U_neighbor[4];
  double numerator, denominator, d2[4], d4[4], area;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  fvm_IO *IO;
  fvm_State FaceState;
  
  Mesh = All->Mesh;
  Data = All->Data;
  IO = All->IO;
  
  fvm_Alloc((void **)&FaceState.W, Mesh->Dim+2, sizeof(double));
  fvm_Alloc((void **)&FaceState.U, Mesh->Dim+2, sizeof(double));
  
  //**************************
  k2 = 1.0/4.0;
  k4 = 3.0/256.0;
  //**************************
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      Data->ElemGroup[egrp].A[elem] = 0.0;
      for (n = 0; n < Mesh->Dim+2; n++)
        Data->ElemGroup[egrp].Nabla2_Q[elem][n] = 0.0;
      numerator =  denominator = 0.0;
      for (f = 0; f < Mesh->ElemGroup[egrp].nFace; f++){
        fgrp = Mesh->ElemGroup[egrp].Face[elem][f].Group;
        face = Mesh->ElemGroup[egrp].Face[elem][f].Number;
        fvm_NeighborAcrossFace(Mesh, egrp, elem, f, &neighbor_grp, &neighbor, &neighbor_face, &flag);
        if (flag == 1){
          area = Data->BFaceGroup[fgrp].Area[face];
          for (k = 0; k < Mesh->Dim; k++)
            Normal[k] = Data->BFaceGroup[fgrp].Normal[face][k];
          for (k = 0; k < Mesh->Dim+2; k++){
            FaceState.W[k] = 0.5*(Data->ElemGroup[egrp].State[elem].W[k]+
                                  Data->BFaceGroup[fgrp].State[face].W[k]);
            Data->ElemGroup[egrp].Nabla2_Q[elem][k] += Data->BFaceGroup[fgrp].State[face].U[k]-
                                                       Data->ElemGroup[egrp].State[elem].U[k];
          }
          
          numerator += fabs(Data->BFaceGroup[fgrp].State[face].W[Mesh->Dim+1]-
                            Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+1]);
          denominator += Data->BFaceGroup[fgrp].State[face].W[Mesh->Dim+1]+
                         Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+1];
        }
        else {
          area = Data->IFace.Area[face];
          for (k = 0; k < Mesh->Dim; k++)
            Normal[k] = Data->IFace.Normal[face][k];
          for (k = 0; k < Mesh->Dim+2; k++){
            FaceState.W[k] = 0.5*(Data->ElemGroup[egrp].State[elem].W[k]+
                                  Data->ElemGroup[neighbor_grp].State[neighbor].W[k]);
            Data->ElemGroup[egrp].Nabla2_Q[elem][k] += Data->ElemGroup[neighbor_grp].State[neighbor].U[k]-
                                                       Data->ElemGroup[egrp].State[elem].U[k];
          }
          
          numerator += fabs(Data->ElemGroup[neighbor_grp].State[neighbor].W[Mesh->Dim+1]-
                            Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+1]);
          denominator += Data->ElemGroup[neighbor_grp].State[neighbor].W[Mesh->Dim+1]+
                         Data->ElemGroup[egrp].State[elem].W[Mesh->Dim+1];
        }
        fvm_ConvertState(FaceState, PRIMITIVE, Mesh->Dim, IO->nChemSpecies);
        Data->ElemGroup[egrp].A[elem] += fabs((FaceState.W[1]*Normal[0] + FaceState.W[2]*Normal[1])*area)+
                                         fvm_SSpeed(FaceState)*area;

      }
      Data->ElemGroup[egrp].nu[elem] = numerator/denominator;
    }
  }
  /****************************************************************************/
  //these two loops have to be separated
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      for (k = 0; k < Mesh->Dim+2; k++) 
        d2[k] = d4[k] = 0.0;
      
      for (f = 0; f < Mesh->ElemGroup[egrp].nFace; f++){
        fgrp = Mesh->ElemGroup[egrp].Face[elem][f].Group;
        face = Mesh->ElemGroup[egrp].Face[elem][f].Number;
        fvm_NeighborAcrossFace(Mesh, egrp, elem, f, &neighbor_grp, &neighbor, &neighbor_face, &flag);
        if (flag == 1){
          nu_neighbor = Data->ElemGroup[egrp].nu[elem];
          A_neighbor = Data->ElemGroup[egrp].A[elem];
          for (k = 0; k < Mesh->Dim+2; k++){
            Nabla2_Q_neighbor[k] = Data->ElemGroup[egrp].Nabla2_Q[elem][k];
            U_neighbor[k] = Data->BFaceGroup[fgrp].State[face].U[k];
          }
        }
        else {
          nu_neighbor = Data->ElemGroup[neighbor_grp].nu[neighbor];
          A_neighbor = Data->ElemGroup[neighbor_grp].A[neighbor];
          for (k = 0; k < Mesh->Dim+2; k++){
            Nabla2_Q_neighbor[k] = Data->ElemGroup[neighbor_grp].Nabla2_Q[neighbor][k];
            U_neighbor[k] = Data->ElemGroup[neighbor_grp].State[neighbor].U[k];
          }
        }
      }
      //Evaluating Epsilon Second
      if (Data->ElemGroup[egrp].nu[elem] > nu_neighbor)
        E2 = k2*Data->ElemGroup[egrp].nu[elem];
      else
        E2 = k2*nu_neighbor;
      //Evaluating Epsilon Fourth
      if ((k4 - E2) > 0.0)
        E4 = k4 - E2;
      else
        E4 = 0.0;
      
      //Evaluating d2 and d4
      for (k = 0; k < Mesh->Dim+2; k++){
        d2[k] += 0.5*(Data->ElemGroup[egrp].A[elem] + A_neighbor)*E2*
                 (U_neighbor[k] - Data->ElemGroup[egrp].State[elem].U[k]);
        d4[k] += 0.5*(Data->ElemGroup[egrp].A[elem] + A_neighbor)*E4*
                 (Nabla2_Q_neighbor[k] - Data->ElemGroup[egrp].Nabla2_Q[elem][k]);
      }
      for (k = 0; k < Mesh->Dim+2; k++)
        Data->ElemGroup[egrp].D[elem][k] = d2[k] - d4[k];
    }
  }
  return 0;
}

/****************************************************************************/
//Function fvm_MFlux
double fvm_MFlux(double p, double rhoQ, double aQ, double pQ)
{
  //This function computes the mass flux between pre and post states.
  double mdot;
  
  //compression
  if(fabs(p/pQ) >= 1.0 - EPS){
    mdot = rhoQ*aQ*sqrt(1.0+(p/pQ - 1.0)*(GAMMA+1.0)/(2.0*GAMMA));
  }
  //expansion
  else{
    mdot = rhoQ*aQ*((GAMMA-1.0)/(2.0*GAMMA))*(1.0-p/pQ)/(1.0-pow((p/pQ),((GAMMA-1.0)/(2*GAMMA))));
  }
  
  return mdot;
}

/****************************************************************************/
//Function fvm_SecRHS
double fvm_SecRHS(double p, double mR, double mL, double uR, double uL, double pR, double pL)
{
  //This function computes the RHS for the root finder method
  double f;
  
  f = p - (mL*pR + mR*pL - mL*mR*(uR-uL))/(mL + mR);
  
  return f;
}

/****************************************************************************/
//Function fvm_Dens
double fvm_Dens(double p, double rhoQ, double pQ)
{
  //This function computes the density for the post states.
  double rho;
  
  //compression
  if(fabs(p/pQ) >= 1.0){
    rho = rhoQ*(1.0+(p/pQ)*(GAMMA+1.0)/(GAMMA-1.0))/((GAMMA+1.0)/(GAMMA-1.0)+(p/pQ));
  }
  //expansion
  else{
    rho = rhoQ*pow((p/pQ),(1.0/GAMMA));
  }
  
  return rho;
}

/****************************************************************************/
//Function fvm_ExpFan
int fvm_ExpFan(double uQ, double aQ, double pQ, double rhoQ, double us, 
               double aQs, double stateI[3], int dir)
{
  //this function computes the primitive variables inside 
  //expansion fan for x/t=0
  double u, rho, p, a;
  if(dir == -1){//left running
    u = uQ - (uQ - aQ)*(us - uQ)/((us - aQs)-(uQ - aQ));
    a = aQ - (uQ - aQ)*(aQs - aQ)/((us - aQs)-(uQ - aQ));
    printf("left: u=%f\n",u);
  }
  else if(dir == 1){//right running
    u = uQ - (uQ + aQ)*(uQ - us)/((uQ + aQ)-(us + aQs));
    a = aQ - (uQ + aQ)*(aQ - aQs)/((uQ + aQ)-(us + aQs));
    printf("right: u=%f\n",u);
  }
  
  p = pQ*pow(a/aQ,2.0/(GAMMA-1.0));
  rho = rhoQ*pow(a/aQ, 2.0*GAMMA/(GAMMA-1.0));
  stateI[0] = rho;
  stateI[1] = u;
  stateI[2] = p;
  
  return 0;
}