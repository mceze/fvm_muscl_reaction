/*
 *  main.c
 *  fvm
 *
 *  Created by Marco Ceze on 7/22/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "main.h"
#include "structure_all.h"
#include "tools_memory.h"
#include "tools_IO.h"
#include "tools_solver.h"


int main()
{
  fvm_All *All;
    
  fvm_CreateAll(&All);
  
  fvm_ReadInput(All);
  
  fvm_Init(All);
  
  fvm_TimeMarching(All);

  fvm_DestroyAll(All);
  
  return 0;
}