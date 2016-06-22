/*
 *  tools_IO.c
 *  fvm
 *
 *  Created by Marco Ceze and Chih-Kuang Kuan on 12/2/09.
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

/****************************************************************************/
//Function ReadInput
int fvm_ReadInput(fvm_All  *All)
{
  FILE *file;
  fvm_IO *IO;
  fvm_Mesh *Mesh;
  int i, j, chem, RefState;
  char string[LONGLEN], dump[LONGLEN], check;
  char temp[LONGLEN];
  //double *P_d, *T_d, *u_d, *v_d, *MolWeight_d, Ta_f_d, Ta_b_d, *HeatF_d, Mref;
  //double Tref, aref, rhoref, *M_d;
  
  IO = All->IO;
  Mesh = All->Mesh;
  
  if((file = fopen("fvm.inp","r")) == NULL){
    printf("No fvm.inp file was found! Exiting!\n");
    exit(0);
  }
  
  /****************************************************************************/
  //Mesh file name
  rewind(file);
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, "MeshFile", 8) == 0) break;
  }
  sscanf(string, "%s %c %s\n", &dump, &check, &(IO->MeshFile));
  
  fvm_ReadMesh(Mesh, IO->MeshFile);
  
  /****************************************************************************/
  //Number of initial states
  rewind(file);
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, "nInitStates", 11) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &(IO->nInitStates));
  
  /****************************************************************************/
  //Reference state
  rewind(file);
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, "RefState", 8) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &RefState);
  RefState--;
  /****************************************************************************/
  //Solve Chemistry?
  rewind(file);
  sprintf(temp, "Chemistry");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %s\n", &dump, &check, &temp);
  if (strncmp(temp, "True", 4) == 0)
    chem = 1;
  else
    chem = 0;
  
  /****************************************************************************/
  //nChemSpecies
  rewind(file);
  sprintf(temp, "nChemSpecies");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &(IO->nChemSpecies));
  
  /****************************************************************************/
  //MolWeight
  rewind(file);
  sprintf(temp, "Molar Weight");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fvm_Alloc((void **)&IO->MolWeight, IO->nChemSpecies, sizeof(double));
  //fvm_Alloc((void **)&MolWeight_d, IO->nChemSpecies, sizeof(double));
  for (i = 0; i < IO->nChemSpecies; i++){
    fscanf(file," %lf",&IO->MolWeight[i]); //Note: This is dimensional!!
  }
  
  /****************************************************************************/
  //HeatOfFormation
  rewind(file);
  sprintf(temp, "Heat Of Formation");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fvm_Alloc((void **)&IO->HeatOfFormation, IO->nChemSpecies, sizeof(double));
  //fvm_Alloc((void **)&HeatF_d, IO->nChemSpecies, sizeof(double));
  for (i = 0; i < IO->nChemSpecies; i++){
    fscanf(file," %lf",&IO->HeatOfFormation[i]);
  }
  
  /****************************************************************************/
  //Arrhenius Parameters
  rewind(file);
  sprintf(temp, "Arrhenius Parameters");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fscanf(file,"%lf %lf %lf\n",&IO->A_f, &IO->b_f, &IO->Ta_f);
  fscanf(file,"%lf %lf %lf\n",&IO->A_b, &IO->b_b, &IO->Ta_b);
  

  /****************************************************************************/
  //Stoichiometric Coefficients
  rewind(file);
  sprintf(temp, "Stoichimetric Coefficients");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fvm_Alloc((void **)&IO->StoichCoeff, IO->nChemSpecies, sizeof(double));
  for (i = 0; i < IO->nChemSpecies; i++){
    fscanf(file," %lf",&IO->StoichCoeff[i]);
  }
  
  /****************************************************************************/
  //Getting initial states
  fvm_Alloc((void **)&IO->InitType, IO->nInitStates, sizeof(int));
  fvm_Alloc2((void ***)&IO->InitGroups, IO->nInitStates, Mesh->nElemGroup, sizeof(int));
  fvm_Alloc((void **)&IO->InitState, IO->nInitStates, sizeof(fvm_State));
  fvm_Alloc((void **)&IO->nEGroups, IO->nInitStates, sizeof(int));
  //dimensional state
  //fvm_Alloc((void **)&P_d, IO->nInitStates, sizeof(double));
  //fvm_Alloc((void **)&T_d, IO->nInitStates, sizeof(double));
  //fvm_Alloc((void **)&u_d, IO->nInitStates, sizeof(double));
  //fvm_Alloc((void **)&v_d, IO->nInitStates, sizeof(double));
  //fvm_Alloc((void **)&M_d, IO->nInitStates, sizeof(double));
  
  for (i = 0; i < IO->nInitStates; i++){
    sprintf(temp, "InitState Block %d", i+1);
    rewind(file);
    while(feof(file) != 1){
      fgets(string, LONGLEN, file);
      if (strncmp(string, temp, strlen(temp)) == 0) break;
    }
    fscanf(file, "%s %c %s\n", &dump, &dump[0], &temp);
    fvm_Alloc((void **)&IO->InitState[i].U, 2 + Mesh->Dim+IO->nChemSpecies, sizeof(double));
    fvm_Alloc((void **)&IO->InitState[i].W, 2 + Mesh->Dim+IO->nChemSpecies, sizeof(double));
    
    if (strcmp(temp, "Conserved") == 0){
      IO->InitType[i] = CONSERVED;
      for (j = 0; j < 2 + Mesh->Dim+IO->nChemSpecies; j++)
        fscanf(file, "%lf", &IO->InitState[i].U[j]);
    }
    else if (strcmp(temp, "Primitive") == 0){
      IO->InitType[i] = PRIMITIVE;
      //fscanf(file, "%lf", &P_d[i]);//dimensional variables
      //fscanf(file, "%lf", &T_d[i]);
      //fscanf(file, "%lf", &u_d[i]);
      //fscanf(file, "%lf", &v_d[i]);
      for (j = 0; j < Mesh->Dim+2+IO->nChemSpecies; j++)
        fscanf(file, "%lf", &IO->InitState[i].W[j]);
    }
    else {
      printf("InitType not reconized! Exiting!\n");
      exit(0);
    }
        
    fscanf(file, "%s %c %d\n",&dump, &dump[0], &IO->nEGroups[i]);
    fscanf(file, "%s %c", &dump, &dump[0]);
    for (j = 0; j < IO->nEGroups[i]; j++){
      fscanf(file, "%d",&IO->InitGroups[i][j]);
      IO->InitGroups[i][j] = IO->InitGroups[i][j]-1;
    }    
  }
  
  //Dimensional average molecular weight for each mixture
  /* for (i = 0; i < IO->nInitStates; i++)
      for (j = 0; j < IO->nChemSpecies; j++)
        M_d[i] = IO->InitState[i].W[4+j]*MolWeight_d[j];
     */
  //Reference average molecular weight
  /* Mref = M_d[RefState];
    
    Tref = T_d[RefState];
    rhoref = P_d[RefState]/((RU/Mref)*Tref);
    aref = sqrt(GAMMA*(RU/Mref)*Tref); */
  
  
  //Non-dimensionalizing the initial states (primitive variables)
  /* for (i = 0; i < IO->nInitStates; i++){
      IO->InitState[i].W[0] = (1.0/rhoref)*P_d[i]/((RU/M_d[i])*T_d[i]);
      IO->InitState[i].W[1] = (1.0/aref)*u_d[i];
      IO->InitState[i].W[2] = (1.0/aref)*v_d[i];
      IO->InitState[i].W[3] = (1.0/(rhoref*pow(aref,2.0)))*P_d[i];
    } */
  
  //Non-dimensionalizing the chemical properties
  /* for (j = 0; j < IO->nChemSpecies; j++){
      IO->MolWeight[j] = MolWeight_d[j]/Mref;
      IO->HeatOfFormation[j] = (1.0/(rhoref*pow(aref,2.0)))*HeatF_d[j];
    }
    IO->Ta_f = Ta_f_d/Tref;
    IO->Ta_b = Ta_b_d/Tref; */
  
  
  /****************************************************************************/
  //CFL
  rewind(file);
  sprintf(temp, "CFL");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %lf\n", &dump, &check, &(IO->CFL));
  
  /****************************************************************************/
  //CalcTime
  rewind(file);
  sprintf(temp, "CalcTime");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %lf\n", &dump, &check, &(IO->CalcTime));
  
  /****************************************************************************/
  //itmax
  rewind(file);
  sprintf(temp, "itmax");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &(IO->itmax));
  
  /****************************************************************************/
  //ConstDt
  rewind(file);
  sprintf(temp, "ConstDt");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &(IO->ConstDt));
  
  /****************************************************************************/
  //Limiter
  rewind(file);
  sprintf(temp, "Limiter");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %s\n", &dump, &check, &temp);
  if (strncmp(temp, "Koren", 5) == 0)
    IO->Limiter = KOREN;
  else if (strncmp(temp, "VanAlbada", 9) == 0)
    IO->Limiter = VANALBADA;
  else if (strncmp(temp, "Superbee", 8) == 0)
    IO->Limiter = SUPERBEE;
  else if (strncmp(temp, "MinMod", 6) == 0)
    IO->Limiter = MINMOD;
  else if (strncmp(temp, "VanLeer", 7) == 0)
    IO->Limiter = VANLEER;
  else if (strncmp(temp, "Ceze", 7) == 0)
    IO->Limiter = CEZE;
  else if (strncmp(temp, "None", 4) == 0)
    IO->Limiter = 0;
  else {
    printf("Limiter not specified. Solver will run in first order mode.\n");
    IO->Limiter = 0;
  }
  
  /****************************************************************************/
  //TimeMarching
  rewind(file);
  sprintf(temp, "TimeMarching");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %s\n", &dump, &check, &temp);
  if (strncmp(temp, "RK3", 5) == 0)
    IO->TimeMarching = RK3;
  else if (strncmp(temp, "FE", 4) == 0)
    IO->TimeMarching = FE;
  else {
    printf("TimeMarching not specified. Solver will run in Forward-Euler mode.\n");
    IO->TimeMarching = FE;
  }
  
  /****************************************************************************/
  //WriteEvery
  rewind(file);
  sprintf(temp, "WriteEvery");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %d\n", &dump, &check, &(IO->WriteEvery));
  
  /****************************************************************************/
  //VelocityReference
  IO->VRef = 1.0; //default value
  rewind(file);
  sprintf(temp, "VelocityReference");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0){
      sscanf(string, "%s %c %lf\n", &dump, &check, &(IO->VRef));
      break; 
    }    
  }
  
  /****************************************************************************/
  //LengthReference
  IO->LRef = 1.0; //default value
  rewind(file);
  sprintf(temp, "LengthReference");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0){
      sscanf(string, "%s %c %lf\n", &dump, &check, &(IO->LRef));
      break; 
    }
  }
  
  /****************************************************************************/
  //Kappa
  IO->Kappa = 1.0/3.0; //default value
  rewind(file);
  sprintf(temp, "Kappa");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0){
      sscanf(string, "%s %c %lf\n", &dump, &check, &(IO->Kappa));
      break;
    }
  }
  
  /****************************************************************************/
  //FluxType
  rewind(file);
  sprintf(temp, "FluxType");
  while(feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  sscanf(string, "%s %c %s\n", &dump, &check, &temp);
  if (strncmp(temp, "RoeE", 4) == 0)
    IO->FluxType = ROEE;
  else if (strncmp(temp, "Jameson", 7) == 0){
    IO->FluxType = JAMESON;
    printf("Time-marching will be set to a 5-stage Runge-Kutta.\n");
  }
  else if (strncmp(temp, "Roe", 3) == 0)
    IO->FluxType = ROE;
  else if (strncmp(temp, "ER", 2) == 0)
    IO->FluxType = ER;
  else {
    printf("FluxType not specified. Solver will run with RoeE.\n");
    IO->FluxType = ROEE;
  }
  /****************************************************************************/
  
  /* fvm_Release((void *)P_d);
    fvm_Release((void *)T_d);
    fvm_Release((void *)u_d);
    fvm_Release((void *)v_d);
    fvm_Release((void *)MolWeight_d); */
  
  return 0;
}

/****************************************************************************/
//Function fvm_ReadMesh
int fvm_ReadMesh(fvm_Mesh *Mesh, char MeshFile[])
{
  FILE *file;
  char string[LONGLEN], dump[LONGLEN], temp[LONGLEN];
  int i, j, k, n, egrp, elem, face, nElemTot, waste, **cells, **eConv;
  int nodes[2], iface, **IntFace;
  file = fopen(MeshFile, "r");
  
  //Reading general mesh information
  rewind(file);
  sprintf(temp, "     NUMNP");
  while (feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  
  fscanf(file, "%d %d %d %d %d %d\n", &Mesh->nNode, &nElemTot, 
         &Mesh->nElemGroup, &Mesh->nBFaceGroup, &Mesh->Dim, &waste);
  
  //Reading the nodes
  rewind(file);
  sprintf(temp, "   NODAL");
  while (feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fvm_Alloc2((void ***)&Mesh->Coord, Mesh->nNode, Mesh->Dim, sizeof(double));
  if (Mesh->Dim == 2)
    for (i = 0; i < Mesh->nNode; i++)
      fscanf(file, "%d %lf %lf\n", &j, &Mesh->Coord[i][0], &Mesh->Coord[i][1]);
  else{
    printf("Dimension not suported!\n");
    exit(0);
  }
  //Reading the cells. OBS: Let's assume quads for now!
  rewind(file);
  sprintf(temp, "      ELEMENTS");
  while (feof(file) != 1){
    fgets(string, LONGLEN, file);
    if (strncmp(string, temp, strlen(temp)) == 0) break;
  }
  fvm_Alloc2((void ***)&cells, nElemTot, 4, sizeof(int));
  for (i = 0; i < nElemTot; i++){
    fscanf(file, "%d %d %d %d %d %d %d\n", &j, &waste, 
           &waste, &cells[i][0], &cells[i][1], 
           &cells[i][2], &cells[i][3]);
    cells[i][0] = cells[i][0]-1;
    cells[i][1] = cells[i][1]-1;
    cells[i][2] = cells[i][2]-1;
    cells[i][3] = cells[i][3]-1;
  }
  //Reading element groups
  fvm_Alloc2((void ***)&eConv, nElemTot, 2, sizeof(int));
  fvm_Alloc((void **)&Mesh->ElemGroup, Mesh->nElemGroup, sizeof(fvm_ElemGroup));
  rewind(file);
  sprintf(temp, "GROUP:");
  for (i = 0; i < Mesh->nElemGroup; i++){
    while (feof(file) != 1){
      fgets(string, LONGLEN, file);
      if (strncmp(string, temp, strlen(temp)) == 0) break;
    }
    sscanf(string, "%s %d %s %d %s %c %s %d", &dump, &j, 
           &dump, &k, &dump, &dump[0], &dump, &n);
    egrp = j-1;
    Mesh->ElemGroup[egrp].nElem = k;
    Mesh->ElemGroup[egrp].nNode = 4;
    Mesh->ElemGroup[egrp].nFace = 4;
    fvm_Alloc2((void ***)&Mesh->ElemGroup[egrp].Node, 
               Mesh->ElemGroup[egrp].nElem, 
               Mesh->ElemGroup[egrp].nNode, sizeof(int));
    fvm_Alloc2((void ***)&Mesh->ElemGroup[egrp].Face, 
               Mesh->ElemGroup[egrp].nElem, 
               Mesh->ElemGroup[egrp].nFace, sizeof(fvm_Face));
    
    fscanf(file, "%s\n", &Mesh->ElemGroup[egrp].Title);
    for (j = 0; j < n; j++) fgets(string, LONGLEN, file);
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      fscanf(file, "%d",&j);
      for (k = 0; k < Mesh->ElemGroup[egrp].nNode; k++)
        Mesh->ElemGroup[egrp].Node[elem][k] = cells[j-1][k];
      //element conversion table
      eConv[j-1][0] = egrp;
      eConv[j-1][1] = elem;
    }
  }
  
  //initializing faces in elements
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      for (face = 0; face < Mesh->ElemGroup[egrp].nFace; face++){
        Mesh->ElemGroup[egrp].Face[elem][face].Number = -1;
        Mesh->ElemGroup[egrp].Face[elem][face].Group = -1;
      }
    }
  }
  
  //Reading boundary faces
  rewind(file);
  fvm_Alloc((void **)&Mesh->BFaceGroup, Mesh->nBFaceGroup, 
            sizeof(fvm_BFaceGroup));
  for (i = 0; i < Mesh->nBFaceGroup; i++){
    sprintf(temp, " BOUNDARY CONDITIONS");
    while (feof(file) != 1){
      fgets(string, LONGLEN, file);
      if (strncmp(string, temp, strlen(temp)) == 0) break;
    }
    fscanf(file, "%s %d %d %d %d\n", &Mesh->BFaceGroup[i].Title, 
           &waste, &Mesh->BFaceGroup[i].nFace, &waste, &waste);
    //filling up the bc type
    if (strcmp(Mesh->BFaceGroup[i].Title, "FARFIELD") == 0)
      Mesh->BFaceGroup[i].Type = FARFIELD;
    if (strcmp(Mesh->BFaceGroup[i].Title, "WALL") == 0)
      Mesh->BFaceGroup[i].Type = WALL;
    if (strcmp(Mesh->BFaceGroup[i].Title, "INFLOWOUTFLOW") == 0)
      Mesh->BFaceGroup[i].Type = INFLOWOUTFLOW;
    if (strcmp(Mesh->BFaceGroup[i].Title, "IN") == 0)
      Mesh->BFaceGroup[i].Type = INFLOWOUTFLOW;
    if (strcmp(Mesh->BFaceGroup[i].Title, "OUT") == 0)
      Mesh->BFaceGroup[i].Type = INFLOWOUTFLOW;
    //filling up the faces
    fvm_Alloc((void **)&Mesh->BFaceGroup[i].BFace, 
              Mesh->BFaceGroup[i].nFace, sizeof(fvm_BFace));
    for (j = 0; j < Mesh->BFaceGroup[i].nFace; j++){
      fscanf(file, "%d %d %d\n", &k, &waste, &n);
      egrp = eConv[k-1][0];
      elem = eConv[k-1][1];
      face = n-1;
      Mesh->BFaceGroup[i].BFace[j].ElemGroup = egrp;
      Mesh->BFaceGroup[i].BFace[j].Elem = elem;
      Mesh->BFaceGroup[i].BFace[j].Face = face;
      Mesh->BFaceGroup[i].BFace[j].Orient = 1;
      Mesh->ElemGroup[egrp].Face[elem][face].Number = j;
      Mesh->ElemGroup[egrp].Face[elem][face].Group = i;
    }
  }
  
  //Internal faces
  fvm_Alloc2((void ***)&IntFace, nElemTot*4, 6, sizeof(int));//this has to change in the future
  
  iface = 0;
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
      for (face = 0; face < Mesh->ElemGroup[egrp].nFace; face++){
        if (Mesh->ElemGroup[egrp].Face[elem][face].Number < 0){
          if (face == Mesh->ElemGroup[egrp].nFace - 1){
            nodes[0] = Mesh->ElemGroup[egrp].Node[elem][Mesh->ElemGroup[egrp].nNode-1];
            nodes[1] = Mesh->ElemGroup[egrp].Node[elem][0];
          }
          else{
            nodes[0] = Mesh->ElemGroup[egrp].Node[elem][face];
            nodes[1] = Mesh->ElemGroup[egrp].Node[elem][face + 1];
          }
          for (i = 0; i < nElemTot; i ++){
            if (nodes[0] == cells[i][1] && nodes[1] == cells[i][0]){
              IntFace[iface][0] = egrp; //ElemGroupL
              IntFace[iface][1] = elem; //ElemL
              IntFace[iface][2] = face; //FaceL
              IntFace[iface][3] = eConv[i][0]; //ElemGroupR
              IntFace[iface][4] = eConv[i][1]; //ElemR
              IntFace[iface][5] = 0; //FaceR
              Mesh->ElemGroup[egrp].Face[elem][face].Number = iface;
              Mesh->ElemGroup[IntFace[iface][3]].Face[IntFace[iface][4]][IntFace[iface][5]].Number = iface;
              iface++;
              break;
            }
            if (nodes[0] == cells[i][2] && nodes[1] == cells[i][1]){
              IntFace[iface][0] = egrp; //ElemGroupL
              IntFace[iface][1] = elem; //ElemL
              IntFace[iface][2] = face; //FaceL
              IntFace[iface][3] = eConv[i][0]; //ElemGroupR
              IntFace[iface][4] = eConv[i][1]; //ElemR
              IntFace[iface][5] = 1; //FaceR
              Mesh->ElemGroup[egrp].Face[elem][face].Number = iface;
              Mesh->ElemGroup[IntFace[iface][3]].Face[IntFace[iface][4]][IntFace[iface][5]].Number = iface;
              iface++;
              break;
            }
            if (nodes[0] == cells[i][3] && nodes[1] == cells[i][2]){
              IntFace[iface][0] = egrp; //ElemGroupL
              IntFace[iface][1] = elem; //ElemL
              IntFace[iface][2] = face; //FaceL
              IntFace[iface][3] = eConv[i][0]; //ElemGroupR
              IntFace[iface][4] = eConv[i][1]; //ElemR
              IntFace[iface][5] = 2; //FaceR
              Mesh->ElemGroup[egrp].Face[elem][face].Number = iface;
              Mesh->ElemGroup[IntFace[iface][3]].Face[IntFace[iface][4]][IntFace[iface][5]].Number = iface;
              iface++;
              break;
            }
            if (nodes[0] == cells[i][0] && nodes[1] == cells[i][3]){
              IntFace[iface][0] = egrp; //ElemGroupL
              IntFace[iface][1] = elem; //ElemL
              IntFace[iface][2] = face; //FaceL
              IntFace[iface][3] = eConv[i][0]; //ElemGroupR
              IntFace[iface][4] = eConv[i][1]; //ElemR
              IntFace[iface][5] = 3; //FaceR
              Mesh->ElemGroup[egrp].Face[elem][face].Number = iface;
              Mesh->ElemGroup[IntFace[iface][3]].Face[IntFace[iface][4]][IntFace[iface][5]].Number = iface;
              iface++;
              break;
            }
            
          }
        }
      }
    }
  }
  
  Mesh->nIFace = iface;
  
  fvm_Alloc((void **)&Mesh->IFace, Mesh->nIFace, sizeof(fvm_IFace));
  
  for(i = 0; i < Mesh->nIFace; i++){
    Mesh->IFace[i].ElemGroupL = IntFace[i][0];
    Mesh->IFace[i].ElemL = IntFace[i][1];
    Mesh->IFace[i].FaceL = IntFace[i][2];
    Mesh->IFace[i].OrientL = 1;
    Mesh->IFace[i].ElemGroupR = IntFace[i][3];
    Mesh->IFace[i].ElemR = IntFace[i][4];
    Mesh->IFace[i].FaceR = IntFace[i][5];
    Mesh->IFace[i].OrientR = 0;
  }
  
  fclose(file);
  fvm_Release2((void **)cells);
  fvm_Release2((void **)eConv);
  fvm_Release2((void **)IntFace);
  
  return 0;
}

/****************************************************************************/
//Function fvm_WriteSolution
int fvm_WriteSolution(fvm_All *All, int it, double time)
{
  int egrp, elem, i, length, nElem, newline_flag, k, j, fgrp, face;
  char outfile[LONGLEN], temp[LONGLEN];
  double TRef;
  FILE *file, *output;
  fvm_IO *IO;
  fvm_Mesh *Mesh;
  fvm_Data *Data;
  
  IO = All->IO;
  Mesh = All->Mesh;
  Data = All->Data;
  
  TRef = IO->LRef/IO->VRef;
  
  time = time/TRef;
  
  length = strlen(IO->MeshFile);
  
  for (i = 0; i < length-4; i++)
    outfile[i] = IO->MeshFile[i];
  
  nElem = 0;
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
    nElem += Mesh->ElemGroup[egrp].nElem;
  
  /****************************************************************************/
  //Writing Tecplot File
  outfile[length-4] = (char)'.';
  outfile[length-3] = (char)'p';
  outfile[length-2] = (char)'l';
  outfile[length-1] = (char)'t';
  outfile[length] = (char)'\0';
  
  if (it == 0){
    file = fopen(outfile, "w");
    output = fopen("rake.m", "w");
  }
  else{
    file = fopen(outfile, "a");
    output = fopen("rake.m", "a");
  }
  
  if (it == 0)
    fprintf(file, "TITLE=\"FVM Solution\"\n");
  
  fprintf(file, "Variables = \"X\" \"Y\" \"RHO\" \"U\" \"V\" \"P\" \"Ya\" \"Yb\" \"Yc\" \n");
  fprintf(file, "ZONE T=\"%d\",SOLUTIONTIME= %lf, N= %d, E= %d, F=FEBLOCK, ET=QUADRILATERAL, VARLOCATION=([3 4 5 6 7 8 9]=CELLCENTERED)\n",it,time,Mesh->nNode, nElem);
  newline_flag = 0;
  j = 0;
  for (i = 0; i < Mesh->nNode; i++){
    fprintf(file, "%lf ",Mesh->Coord[i][0]);
    j++;
    if (j/8 > newline_flag){
      newline_flag++;
      fprintf(file,"\n");
    }
  }
  
  for (i = 0; i < Mesh->nNode; i++){
    fprintf(file, "%lf ",Mesh->Coord[i][1]);
    j++;
    if (j/8 > newline_flag){
      newline_flag++;
      fprintf(file,"\n");
    }
  }
  for (k = 0; k < Mesh->Dim+2+IO->nChemSpecies; k++){
    for (egrp = 0; egrp < Mesh->nElemGroup; egrp++){
      for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++){
        fprintf(file, "%lf ", Data->ElemGroup[egrp].State[elem].W[k]);
        j++;
        if (j/8 > newline_flag){
          newline_flag++;
          fprintf(file,"\n");
        }
      }
    }
  }
  
  fprintf(file, "\n");
  
  for (egrp = 0; egrp < Mesh->nElemGroup; egrp++)
    for (elem = 0; elem < Mesh->ElemGroup[egrp].nElem; elem++)
      fprintf(file, "%d %d %d %d\n", Mesh->ElemGroup[egrp].Node[elem][0]+1, 
              Mesh->ElemGroup[egrp].Node[elem][1]+1, Mesh->ElemGroup[egrp].Node[elem][2]+1, 
              Mesh->ElemGroup[egrp].Node[elem][3]+1);
  
  
  fclose(file);
  
  /****************************************************************************/
  //Rake file
  
  if (it == 0)
    fprintf(output, "%% Rake at OUT Boundary\nclear all\n");
  
  fprintf(output, "t(%d) = %lf;\n", (it/IO->WriteEvery)+1, time);
  sprintf(temp, "OUT");
  for (fgrp = 0; fgrp < Mesh->nBFaceGroup; fgrp++){
    if (strcmp(Mesh->BFaceGroup[fgrp].Title, temp) == 0){
      for (face = 0; face < Mesh->BFaceGroup[fgrp].nFace; face++){
        elem = Mesh->BFaceGroup[fgrp].BFace[face].Elem;
        egrp = Mesh->BFaceGroup[fgrp].BFace[face].ElemGroup;
        fprintf(output, "output(%d,:,%d) = [", face+1, (it/IO->WriteEvery)+1);
        fprintf(output, "%lf %lf %lf %lf %lf %lf];\n", Data->ElemGroup[egrp].ElemCenter[elem][0],
                Data->ElemGroup[egrp].ElemCenter[elem][1], Data->ElemGroup[egrp].State[elem].W[0],
                Data->ElemGroup[egrp].State[elem].W[1], Data->ElemGroup[egrp].State[elem].W[2], 
                Data->ElemGroup[egrp].State[elem].W[3]);
      }
    }
  }
  
  fclose(output);
  
  return 0;
}



