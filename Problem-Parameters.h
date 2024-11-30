/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#pragma once
//#define PP_DEBUG
//=========================== Method Parameters =========================
#define PP_EPS_ZERO		1E-9	// Accuracy for comparison with zero
#define PP_ETA_TO_APEX	1	// Distance from apex base to apex point (DEFAULT: 1000000)
/*============================== rnd5-1-2 LP problem ==============================*/
//#define PP_PROBLEM_NAME	"rnd5-1-2"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
//------------------------------------------------------------------/**/
#define PP_MM (2*PP_M+2*PP_N)
#define PP_MAX_NUM_SHIFTS_SAME_LENGTH	5 // Maximal number of shifts with the same length
#define PP_MAX_ITER_COUNT				10000000000 // Maximal count of iterations
#define PP_ADD_FLAG						PP_N
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_DATABASE_OUTPUT
#define PP_OUTPUT_LIMIT	30	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 14
//------------------------- Matrix format ----------------
#define PP_INFINITY				1E+308		// Highest bound in *_hi.mtx
//-------------------------- Jobs  -----------------------
#define PP_JOB_PSEUDOPOJECTION	0 
#define PP_JOB_LEASTPOJECTION	2 
//#define PP_USE_LEASTPROJECTION
//-------------------------- Process states --------------------------
#define PP_STATE_START						0
#define PP_STATE_FIND_INTERIOR_POINT		1
#define PP_STATE_FIND_INITIAL_APPROXIMATION	2
//------------- Vector Projection Onto Halfspace exit codes -------------
#define PP_EXITCODE_DEGENERATE_INEQUALITY	0
#define PP_EXITCODE_INSIDE_HALFSPACE		1
#define PP_EXITCODE_SUCCESSFUL_PROJECTING	2
#define PP_EXITCODE_PARALLEL				3
#define PP_EXITCODE_RECESSIVE				4
#define PP_EXITCODE_ON_HYPERPLANE			5

//=========================== Problem Parameters =========================
#define PP_MAX_RND_INEQUALITIES	400
#define PP_MAX_N				400
#define PP_MAX_MTX_N			800
#define PP_MAX_M				400
#define PP_MAX_MTX_M			800

#define PP_FILE_INI "config.ini"
static std::string PP_PATH;
//static int PP_N;												// Space dimension
