/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: BSF-Types.h (Problem Independent Types)
Prefix: BT
Author: Nikolay A. Olkhovsky 
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#pragma once
#include "Problem-bsfTypes.h"	// Predefined BSF Problem Types
//=========================== BSF Types ===============================
struct BT_order_T {
	char exit;		// true, if worker must stop
	int jobCase;
	int iterCounter;
	PT_bsf_parameter_T parameter;
};

struct BT_extendedReduceElem_T {// Extended element type of reduce list
	PT_bsf_reduceElem_T elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_1 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_1 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_2 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_2 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_3 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_3 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};