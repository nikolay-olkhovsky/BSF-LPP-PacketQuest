/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: Problem-bsfTypes.h (Predefined Problem-depended LBSF Types)
Prefix: PT_bsf
Author: Nikolay A. Olkhovsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {
	PT_vector_T x;			// Current approximation
	PT_float_T max_residual;// Maximum of residuals
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	PT_float_T* a;				// Pointer to row: a[0],...,a[n-1]
	PT_float_T* b;				// Pointer to constant term b: a_0x_0+...+a_{n-1}x_{n-1} \leq b
};

// 0. Pseudoprojection
struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T projection;		// Orthogonal projection onto hyperplane
	int nonZeroCounter;			// Counter of nonzero elements
	PT_float_T residual;		// <a,x> - b
	bool activeTag[PP_MM];		// Tag of activity
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

// 2. Least projection
struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	PT_vector_T projection;		// Orthogonal projection onto hyperplane
	int nonZeroCounter;			// Counter of nonzero elements
	PT_float_T residual;		// <a,x> - b
	bool activeTag[PP_MM];		// Tag of activity
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Not used
};