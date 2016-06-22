/*
 *  tools_solver.h
 *  fvm
 *
 *  Created by Marco Ceze on 8/3/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#ifndef _fvm_tools_solver_h
#define _fvm_tools_solver_h 1

#include "structure_mesh.h"
#include "structure_data.h"
#include "structure_IO.h"
#include "structure_all.h"

/****************************************************************************/
//Function fvm_CalculateFaceData
int fvm_CalculateFaceData(fvm_All *All);

/****************************************************************************/
//Function fvm_ConvertState
int fvm_ConvertState(fvm_State State, int Source, int Dim, int  nChemSpecies);

/****************************************************************************/
//Function fvm_Init
int fvm_Init(fvm_All *All);

/****************************************************************************/
//Function fvm_SSpeed
double fvm_SSpeed(fvm_State State);

/****************************************************************************/
//Function fvm_CompFlux
int fvm_CompFlux(fvm_State State, double F[4], double G[4]);

/****************************************************************************/
//Function fvm_FluxRoe (Still 2D)
int fvm_FluxRoe(fvm_All *All, fvm_State StateL, fvm_State StateR, 
                fvm_Face Face, double NormalFlux[7], double *MaxWaveSpeed);

/****************************************************************************/
//Function fvm_CompFaceStates
int fvm_CompFaceStates(fvm_All *All, fvm_Face Face, fvm_State StateL, fvm_State StateR);

/****************************************************************************/
//Function fvm_FindOppositeFace
int fvm_FindOppositeFace(int face, int egrp, fvm_All *All);

/****************************************************************************/
//Function fvm_NeighborAcrossFace
int fvm_NeighborAcrossFace(fvm_Mesh *Mesh, int egrp, int elem, int face, int *neighbor_grp, 
                           int *neighbor, int *neighbor_face, int *is_it_boundary);

/****************************************************************************/
//Funtion fvm_KorenPhi
double fvm_KorenPhi(double r);

/****************************************************************************/
//Function fvm_ApplyBC
int fvm_ApplyBC(fvm_Face Face, fvm_All *All, fvm_State InsideState, fvm_State OutsideState);

/****************************************************************************/
//Function fvm_CompEulerRHS
int fvm_CompEulerRHS(fvm_All *All);

/****************************************************************************/
//Function fvm_CompChemRHS
int fvm_CompChemRHS(fvm_All *All);

/****************************************************************************/
//Function fvm_TimeMarching
int fvm_TimeMarching(fvm_All *All);

/****************************************************************************/
//Function fvm_CalcStepSize
int fvm_CalcStepSize(fvm_All *All);

/****************************************************************************/
//Function fvm_VanAlbadaPhi
double fvm_VanAlbadaPhi(double r);

/****************************************************************************/
//Function fvm_SuperbeePhi
double fvm_SuperbeePhi(double r);

/****************************************************************************/
//Function fvm_VanLeerPhi
double fvm_VanLeerPhi(double r);

/****************************************************************************/
//Function fvm_CezePhi
double fvm_CezePhi(double r);

/****************************************************************************/
//Function fvm_CalculatePhi
double fvm_CalculatePhi(fvm_IO *IO, double r);

/****************************************************************************/
//Function fvm_Convective
int fvm_Convective(fvm_All *All);

/****************************************************************************/
//Function fvm_BCond
int fvm_BCond(fvm_All *All);

/****************************************************************************/
//Function fvm_Dissipative
int fvm_Dissipative(fvm_All *All);

/****************************************************************************/
//Function fvm_MFlux
double fvm_MFlux(double p, double rhoQ, double aQ, double pQ);

/****************************************************************************/
//Function fvm_SecRHS
double fvm_SecRHS(double p, double mR, double mL, double uR, double uL, double pR, double pL);

/****************************************************************************/
//Function fvm_Dens
double fvm_Dens(double p, double rhoQ, double pQ);

/****************************************************************************/
//Function fvm_ExpFan
int fvm_ExpFan(double uQ, double aQ, double pQ, double rhoQ, double us, 
               double aQs, double stateI[3], int dir);

/****************************************************************************/
//Function fvm_FluxER (Still 2D)
int fvm_FluxER(fvm_All *All, fvm_State StateL, fvm_State StateR, 
                fvm_Face Face, double NormalFlux[4], double *MaxWaveSpeed);

#endif