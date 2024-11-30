/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: Problem-Types.h (LBSF Types)
Prefix: PT
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef double PT_float_T;
typedef unsigned long long PT_unsigned_T;
typedef PT_float_T	PT_vector_T[PP_MAX_N];
typedef PT_float_T	PT_matrix_T[PP_MAX_M][PP_MAX_N];
typedef PT_float_T	PT_column_T[PP_MAX_M];
#ifdef PP_DATABASE_OUTPUT
struct Problem {
    unsigned id;
    int N;
    int seed;
    double high;
    double low;
    std::vector<char> c;
};

struct Inequality {
    unsigned id;
    std::vector<char> coefficients;
    double b;
    int problem_id;
};

struct SurfacePoint {
    unsigned id;
    std::vector<char> coefficients;
    int problem_id;
};

struct Precedent {
    unsigned id;
    int problem_id;
    std::vector<char> coefficients;
    int face;
    std::vector<char> d;
    double shift;
    std::vector<char> face_count;
    std::vector<char> face_numbers;
};

struct Image {
    unsigned id;
    int precedent_id;
    double density;
    double rank;
    std::vector<char> answer_vector;
    std::vector<char> cosine_vector;
    int num_of_points;
    std::vector<char> data;
    std::optional<std::vector<char>> field_points;
};
#endif // PP_DATABASE_OUTPUT