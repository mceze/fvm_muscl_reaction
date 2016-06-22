/*
 *  main.h
 *  fvm
 *
 *  Created by Marco Ceze on 7/22/09.
 *  Copyright 2009 University of Michigan. All rights reserved.
 *
 */

#ifndef _fvm_main_h
#define _fvm_main_h 1

/****************************************************************************/
//Definitions
#define GAMMA 1.4
#define RU 8.314
#define EPS 1.0e-8
#define SHORTLEN 30
#define MEDIUMLEN 50
#define LONGLEN 100
#define CONSERVED -1
#define PRIMITIVE -2
#define INTERIORFACE -1
#define RSitmax 100

//Flux Functions Macros
#define ROEE 1
#define ROE 2
#define ER 3
#define JAMESON 4

//Direction macros
#define RIGHT 1
#define LEFT -1

//Limiter Macros
#define KOREN 1
#define VANALBADA 2
#define SUPERBEE 3
#define MINMOD 4
#define VANLEER 5
#define CEZE 6

//TimeMarching Macros
#define RK3 1
#define FE 2

//Boundary Condition types
#define INTERIOR 1
#define WALL 2
#define FARFIELD 3
#define INFLOWOUTFLOW 4

/****************************************************************************/

#endif