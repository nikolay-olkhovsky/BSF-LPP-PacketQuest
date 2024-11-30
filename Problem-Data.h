/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;

static PT_unsigned_T PD_index = 0;			// Index of current LPP in dataset
static PT_unsigned_T PD_packetSize;         // Number of LPPs in dataset
static double        PD_time;               // Overall time of execution

//========================== Problem variables ====================================
static int		PD_m;						// Current number of inequalities
static int		PD_n;						// Space dimension
static int		PD_m_sp;					// Number of hyperplanes that include surface point
static int		PD_numShiftsSameLength;		// Number of shifts with the same length
static bool		PD_pointIn;					// Point is inside polytope
static int		PD_state;					// State of Job Dispatcher (see PC_bsf_JobDispatcher)
static double	PD_objF_u;
//========================== Problem structures ====================================
static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities 
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;					// Objective Function Coefficients
static PT_vector_T PD_apexPoint;			// Apex point
static PT_vector_T PD_x0;					// Starting point (defauil is 0)
static PT_vector_T PD_u;					// Base point on Polytope
static PT_vector_T PD_direction;			// Unit vector to set shift direction
static PT_vector_T PD_hi;					// Higher bound
static PT_vector_T PD_lo;					// Lower bound
static PT_vector_T PD_e_c;					// e_c = c/||c||
static bool PD_active_tag[PP_MM];			// Tag of activity during least projecting
static bool PD_recessive_tag[PP_MM];		// Tag of recessivity
static bool PD_u_hyperplanes_tag[PP_MM];	// Tag of point u belongness to recessive hyperplane

//========================== Files ==============================
#ifdef PP_DATABASE_OUTPUT
static auto storage = sqlite_orm::make_storage("C:/HS/dataset100000_3.sqlite3",
    sqlite_orm::make_table("problems",
        sqlite_orm::make_column("id", &Problem::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("N", &Problem::N),
        sqlite_orm::make_column("seed", &Problem::seed),
        sqlite_orm::make_column("high", &Problem::high),
        sqlite_orm::make_column("low", &Problem::low),
        sqlite_orm::make_column("c", &Problem::c)
    ),
    sqlite_orm::make_table("inequalities",
        sqlite_orm::make_column("id", &Inequality::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &Inequality::coefficients),
        sqlite_orm::make_column("b", &Inequality::b),
        sqlite_orm::make_column("problem_id", &Inequality::problem_id),
        sqlite_orm::foreign_key(&Inequality::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("surface_points",
        sqlite_orm::make_column("id", &SurfacePoint::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &SurfacePoint::coefficients),
        sqlite_orm::make_column("problem_id", &SurfacePoint::problem_id),
        sqlite_orm::foreign_key(&SurfacePoint::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("precedents",
        sqlite_orm::make_column("id", &Precedent::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("coefficients", &Precedent::coefficients),
        sqlite_orm::make_column("face", &Precedent::face),
        sqlite_orm::make_column("d", &Precedent::d),
        sqlite_orm::make_column("shift", &Precedent::shift),
        sqlite_orm::make_column("face_count", &Precedent::face_count),
        sqlite_orm::make_column("face_numbers", &Precedent::face_numbers),
        sqlite_orm::make_column("problem_id", &Precedent::problem_id),
        sqlite_orm::foreign_key(&Precedent::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("images",
        sqlite_orm::make_column("id", &Image::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("density", &Image::density),
        sqlite_orm::make_column("rank", &Image::rank),
        sqlite_orm::make_column("answer_vector", &Image::answer_vector),
        sqlite_orm::make_column("cosine_vector", &Image::cosine_vector),
        sqlite_orm::make_column("num_of_points", &Image::num_of_points),
        sqlite_orm::make_column("data", &Image::data),
        sqlite_orm::make_column("precedent_id", &Image::precedent_id),
        sqlite_orm::make_column("field_points", &Image::field_points),
        sqlite_orm::foreign_key(&Image::precedent_id).references(&Precedent::id)
    )
);
std::vector<unsigned> PD_ids;
std::vector<Problem> PD_problems;
std::vector<Inequality> PD_inequalities;
#else
CMTXX0PacketWriter* PD_packetWriter;
CMTXReader* PD_packetReader;
CProblem* PD_currentProblem;
#endif // PP_DATABASE_OUTPUT
