/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP starting point initializer (Quest)
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
// PP_STATE_START
// PP_STATE_FIND_INITIAL_APPROXIMATION
// PP_STATE_DETERMINE_DIRECTION

#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PD_n; j++) // Generating initial approximation
		parameter->x[j] = PD_u[j];
};

void PC_bsf_Start(bool* success) {
	ini::IniFile config;
	config.load(PP_FILE_INI);
	PP_PATH = config["general"]["PP_PATH"].as<string>();
//	PP_N = config["general"]["PP_N"].as<int>();

#ifdef PP_DATABASE_OUTPUT
	storage.sync_schema(true);
	PD_ids = storage.select(&Problem::id);
	PD_packetSize = PD_ids.size();
	PD_problems = storage.get_all<Problem>();
	PD_inequalities = storage.get_all<Inequality>();

	storage.begin_transaction();
#else

	PD_packetReader = new CMTXReader(PP_PATH.c_str());
	PD_packetSize = PD_packetReader->packetSize;
	if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
	{
		PD_packetWriter = new CMTXX0PacketWriter(PP_PATH.c_str(), PD_packetReader->packetSize);
		//PD_packetWriter->clearFolder();
		PD_packetWriter->open();
	}
#endif // PP_DATABASE_OUTPUT
	PD_index = 0;
	PD_time = clock();
}

void PC_bsf_Init(bool* success) {
	PD_state = PP_STATE_START;
//	PD_problemName = PP_PROBLEM_NAME;

//	*success = LoadMatrixFormat();
//	if (*success == false)
//		return;
#ifdef PP_DATABASE_OUTPUT
	vector<double> _buff;
	vector<Inequality> inequalities;
	//auto problem = storage.get<Problem>(PD_ids[PD_index]);
	auto problem = PD_problems[PD_index];
	//auto inequalities = storage.get_all<Inequality>(sqlite_orm::where(sqlite_orm::c(&Inequality::problem_id) == problem.id));
	copy_if(PD_inequalities.begin(), PD_inequalities.end(), back_inserter(inequalities), [&problem](const Inequality& item)->bool {return item.problem_id == problem.id; });

	PD_n = problem.N;
	PD_m = 0;

	for (unsigned i = 0; i < PD_n; i++) {
		for (unsigned j = 0; j < PD_n; j++)
			if (i == j) PD_A[PD_m][j] = 1.;
			else		PD_A[PD_m][j] = 0.;
		PD_b[PD_m] = problem.high;
		PD_m++;
	}
	for (auto& inequality : inequalities) {
		_buff = charToDouble(inequality.coefficients);
		for (unsigned j = 0; j < PD_n; j++)
			PD_A[PD_m][j] = _buff[j];
		PD_b[PD_m] = inequality.b;
		PD_m++;
	}
	for (unsigned i = 0; i < PD_n; i++) {
		for (unsigned j = 0; j < PD_n; j++)
			if (i == j) PD_A[PD_m][j] = -1.;
			else		PD_A[PD_m][j] = 0.;
		PD_b[PD_m] = problem.low;
		PD_m++;
	}

	_buff = charToDouble(problem.c);
	for (unsigned i = 0; i < PD_n; i++) {
		PD_c[i] = _buff[i];
	}
#else
	PD_currentProblem = PD_packetReader->readProblem();

	PD_m = PD_currentProblem->A->getRows();
	PD_n = PD_currentProblem->A->getCols();

	for (unsigned i = 0; i < PD_m; i++) {
		for (unsigned j = 0; j < PD_n; j++)
			PD_A[i][j] = PD_currentProblem->A->getValue(i, j);
		PD_b[i] = PD_currentProblem->b->getValue(i);
	}
	for (unsigned i = 0; i < PD_n; i++) {
		PD_lo[i] = PD_currentProblem->lo->getValue(i);
		PD_c[i] = -PD_currentProblem->c->getValue(i);
		PD_hi[i] = PD_currentProblem->hi->getValue(i);
	}

	*success = Conversion();
	if (*success == false)
		return;
#endif // PP_DATABASE_OUTPUT

	PD_index++;
	for (int j = 0; j < PD_n; j++){
		PD_u[j] = 0;
		PD_x0[j] = 0;
	}

	UnitObjVector(PD_e_c);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_MM;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);
}

// 0. Pseudo-pojection
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
) {
	// List element no: BSF_sv_addressOffset+BSF_sv_numberInSublist
	int exitCode;

	for (int i = 0; i < PD_m; i++)
		reduceElem->activeTag[i] = false;

	reduceElem->residual =
		Vector_OrthogonalProjectionOntoHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection, &exitCode);

	switch (exitCode) {
	case PP_EXITCODE_DEGENERATE_INEQUALITY:
	case PP_EXITCODE_INSIDE_HALFSPACE:
		reduceElem->nonZeroCounter = 0;
		return;
	case PP_EXITCODE_SUCCESSFUL_PROJECTING:
		if (fabs(reduceElem->residual) < PP_EPS_ZERO && BSF_sv_parameter.max_residual < PP_EPS_ZERO) {
			reduceElem->nonZeroCounter = 0;
			return;
		}
		break;
	default:
		cout << "Vector_ProjectOnHalfspace: Undefined exit code!" << endl;
		*success = false;
		return;
	}
	reduceElem->nonZeroCounter = 1;
	reduceElem->activeTag[BSF_sv_addressOffset + BSF_sv_numberInSublist] = true;
}

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// not used
}

// 2. Least projection
void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// number of hyperplane = BSF_sv_addressOffset + BSF_sv_numberInSublist
	int exitCode;
	PT_vector_T o; // vector of oblique projection
	PT_float_T NormSquare_r;
	PT_float_T NormSquare_o;

	for (int i = 0; i < PD_m; i++)
		reduceElem->activeTag[i] = false;

	//cout << "Recessive tag = " << PD_recessive_tag[BSF_sv_addressOffset + BSF_sv_numberInSublist] << endl;
	if (!PD_recessive_tag[BSF_sv_addressOffset + BSF_sv_numberInSublist]) {
		reduceElem->nonZeroCounter = 0;
		return;
	}

	reduceElem->residual =
		Vector_OrthogonalProjectionOntoHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection, &exitCode);

	cout << "Vector_OrthogonalProjectionOntoHalfspace: " << exitCode << endl;
	switch (exitCode) {
	case PP_EXITCODE_DEGENERATE_INEQUALITY:
	case PP_EXITCODE_INSIDE_HALFSPACE:
		reduceElem->nonZeroCounter = 0;
		return;
	case PP_EXITCODE_SUCCESSFUL_PROJECTING:
		if (fabs(reduceElem->residual) < PP_EPS_ZERO && BSF_sv_parameter.max_residual < PP_EPS_ZERO &&
			BSF_sv_parameter.max_residual != -PP_INFINITY) {
			reduceElem->nonZeroCounter = 0;
			return;
		}
		break;
	default:
		cout << "Vector_OrthogonalProjectionOntoHalfspace: Undefined exit code!" << endl;
		*success = false;
		return;
	}

	NormSquare_r = Vector_NormSquare(reduceElem->projection);

	for (int i = 0; i < PD_m; i++) {

		if (fabs(Vector_DotProduct(PD_A[i], BSF_sv_parameter.x) - PD_b[i]) < PP_EPS_ZERO && BSF_sv_parameter.max_residual < PP_EPS_ZERO
			&& BSF_sv_parameter.max_residual != -PP_INFINITY)
			continue;

		Vector_ObliqueProjectionOntoHalfspace(BSF_sv_parameter.x, PD_A[i], PD_b[i], reduceElem->projection, o, &exitCode);

		cout << "Vector_ObliqueProjectionOntoHalfspace: " << exitCode << endl;
		switch (exitCode) {
		case PP_EXITCODE_INSIDE_HALFSPACE:
			continue;
		case PP_EXITCODE_ON_HYPERPLANE:
			continue;
		case PP_EXITCODE_RECESSIVE:
			continue;
		case PP_EXITCODE_PARALLEL:
			continue;
		case PP_EXITCODE_SUCCESSFUL_PROJECTING:
			break;
		default:
			cout << "Vector_ObliqueProjectionOntoHalfspace: Invalid exit code!" << endl;
			assert(false);
		}

		NormSquare_o = Vector_NormSquare(o);
		if (NormSquare_o < NormSquare_r - PP_EPS_ZERO) {
			for (int j = 0; j < PD_n; j++)
				reduceElem->projection[j] = 0;
			reduceElem->nonZeroCounter = 0;
			reduceElem->residual = -1;
			reduceElem->activeTag[BSF_sv_addressOffset + BSF_sv_numberInSublist] = false;
			return;
		}
	}

	reduceElem->nonZeroCounter = 1;
	reduceElem->activeTag[BSF_sv_addressOffset + BSF_sv_numberInSublist] = true;
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudoprojection
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	Vector_Addition(x->projection, y->projection, z->projection);
	z->nonZeroCounter = x->nonZeroCounter + y->nonZeroCounter;
	z->residual = PF_MAX(x->residual, y->residual);
	for (int i = 0; i < PD_m; i++) {
		z->activeTag[i] = x->activeTag[i] || y->activeTag[i];
	}
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// not used
}

// 2. Least projection
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	Vector_Addition(x->projection, y->projection, z->projection);
	z->nonZeroCounter = x->nonZeroCounter + y->nonZeroCounter;
	z->residual = PF_MAX(x->residual, y->residual);
	for (int i = 0; i < PD_m; i++) {
		z->activeTag[i] = x->activeTag[i] || y->activeTag[i];
	}
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

// 0. Pseudoprojection
void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->nonZeroCounter == 0;
	if (PD_pointIn)
		return;

	parameter->max_residual = reduceResult->residual;

	for (int i = 0; i < PD_m; i++)
		PD_active_tag[i] = reduceResult->activeTag[i];

	static PT_vector_T relaxationVector;
	for (int j = 0; j < PD_n; j++)
		relaxationVector[j] = reduceResult->projection[j] / (double)(reduceResult->nonZeroCounter);
	Vector_PlusEquals(parameter->x, relaxationVector);
}

// 1. Movement on Polytope  ========================================================
void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults_2: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->nonZeroCounter == 0;
	cout << "Inner point found: " << PD_pointIn << endl;
	if (PD_pointIn)
		return;

	parameter->max_residual = reduceResult->residual;

	for (int i = 0; i < PD_m; i++)
		PD_active_tag[i] = reduceResult->activeTag[i];

	static PT_vector_T relaxationVector;
	for (int j = 0; j < PD_n; j++)
		relaxationVector[j] = reduceResult->projection[j] / (double)(reduceResult->nonZeroCounter);
	Vector_PlusEquals(parameter->x, relaxationVector);
}

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit,
	double t
) {
	static PT_vector_T shiftBasePoint;
	static double* ptr_unitVectorToSurface;
//	const char* x0_File = PD_MTX_File_x0.c_str();

	switch (PD_state) {
	case PP_STATE_START://-------------------------- Start -----------------------------
		if (PointInPolytope(PD_x0)) {
			Vector_Copy(PD_x0, PD_apexPoint);
			ApexPoint(PD_x0, PD_apexPoint);
#ifdef PP_DEBUG
			cout << "Apex point:\t";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << PD_apexPoint[j];
			if (PP_OUTPUT_LIMIT < PD_n)
				cout << " ...";
			cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_apexPoint);
			cout << endl;

			for (int i = 0; i < PD_m; i++)
				if (!PD_recessive_tag && PointInHalfspace(PD_apexPoint, PD_A[i], PD_b[i])) {
					if (Vector_OnHyperplane(PD_apexPoint, PD_A[i], PD_b[i]))
						continue;
					cout << "Apex point inside half-space " << i << "!" << endl;
				}

			cout << "Apex point belongs to hyperplane: ";
			for (int i = 0; i < PD_m; i++) {
				if (Vector_OnHyperplane(PD_apexPoint, PD_A[i], PD_b[i]))
					cout << i << " ";
			}
			cout << endl;

			cout << "--------- Finding initial_approximation ------------\n";
#endif
			// Preparations for finding initial approximation
			Vector_Copy(PD_apexPoint, parameter->x);
#ifdef PP_USE_LEASTPROJECTION
			* job = PP_JOB_LEASTPOJECTION;
			parameter->max_residual = -PP_INFINITY;
#else
			* job = PP_JOB_PSEUDOPOJECTION;
#endif
			PD_state = PP_STATE_FIND_INITIAL_APPROXIMATION;
			break;
		}
		// Preparations for finding interior point
		Vector_Copy(PD_x0, parameter->x);
		PD_state = PP_STATE_FIND_INTERIOR_POINT;
#ifdef PP_DEBUG
		cout << "--------- Finding interior point ------------\n";
#endif
		break;
	case PP_STATE_FIND_INTERIOR_POINT://------------------------- Finding interior point -----------------------------
		if (!PD_pointIn) {
			return;
		}
		ApexPoint(parameter->x, PD_apexPoint);
#ifdef PP_DEBUG
		cout << "Apex point:\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << PD_apexPoint[j];
		if (PP_OUTPUT_LIMIT < PD_n)
			cout << " ...";
		cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_apexPoint);
		cout << endl;
		cout << "--------- Finding initial_approximation ------------\n";
#endif
		// Preparations for finding initial approximation
		Vector_Copy(PD_apexPoint, parameter->x);
#ifdef PP_USE_LEASTPROJECTION
		* job = PP_JOB_LEASTPOJECTION;
		parameter->max_residual = -PP_INFINITY;
#else
		* job = PP_JOB_PSEUDOPOJECTION;
#endif 
		PD_state = PP_STATE_FIND_INITIAL_APPROXIMATION;
		break;
	case PP_STATE_FIND_INITIAL_APPROXIMATION://-------------------------- Finding initial approximationt -----------------------------
		if (!PD_pointIn) {
			return;
		}

		Vector_Copy(parameter->x, PD_u);
		PD_objF_u = ObjF(PD_u);

		PD_m_sp = 0;
		for (int i = 0; i < PD_m; i++) {
			if (Vector_OnHyperplane(PD_u, PD_A[i], PD_b[i])) {
				PD_u_hyperplanes_tag[i] = true;
				PD_m_sp++;
			}
			else
				PD_u_hyperplanes_tag[i] = false;
		}

#ifdef PP_DEBUG
		cout << "u belongs to hyperplanes: ";
		for (int i = 0; i < PD_m; i++) {
			if (PD_u_hyperplanes_tag[i])
				cout << i << " ";
		}
		cout << endl;
#endif // PP_DEBUG

//		* exit = true;
		if (PD_index == PD_packetSize) {
#ifdef PP_DATABASE_OUTPUT
			storage.commit();
#endif //PP_DATABASE_OUTPUT
			*exit = true;
			return;
		}
		*job = BD_JOB_RESET;
		break;

	default://------------------------------------- default -----------------------------------
		cout << "PC_bsf_JobDispatcher: Undefined state!" << endl;
//		*exit = true;
		if (PD_index == PD_packetSize) {
#ifdef PP_DATABASE_OUTPUT
			storage.commit();
#endif //PP_DATABASE_OUTPUT
			*exit = true;
			return;
		}
		*job = BD_JOB_RESET;
		break;
	}
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
//	cout << "=================================================== Quest ====================================================" << endl;
////	cout << "Problem name: " << PD_problemName << endl;
//#ifdef PP_USE_LEASTPROJECTION
//	cout << "Mode: Least projections " << endl;
//#else
//	cout << "Mode: Pseudoprojections " << endl;
//#endif // PP_USE_LEASTPROJECTION
//	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
//#ifdef PP_BSF_OMP
//#ifdef PP_BSF_NUM_THREADS
//	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
//#else
//	cout << "Number of Threads: " << omp_get_num_procs() << endl;
//#endif // PP_BSF_NUM_THREADS
//#else
//	cout << "OpenMP is turned off!" << endl;
//#endif // PP_BSF_OMP
//
//#ifdef PP_BSF_FRAGMENTED_MAP_LIST
//	cout << "Map List is Fragmented." << endl;
//#else
//	cout << "Map List is not Fragmented." << endl;
//#endif
//	cout << "Before conversion: m =\t" << PP_M << "\tn = " << PP_N << endl;
//	cout << "After conversion:  m =\t" << PD_m << "\tn = " << PD_n << endl;
//	cout << "Eps Zero: " << PP_EPS_ZERO << endl;
//	cout << "Eta to apex: " << PP_ETA_TO_APEX << endl;
//	cout << "--------------- Data ---------------\n";
//#ifdef PP_MATRIX_OUTPUT
//	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
//	for (int i = 0; i < PD_m; i++) {
//		cout << i << ")";
//		for (int j = 0; j < PD_n; j++)
//			cout << setw(PP_SETW) << PD_A[i][j];
//		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
//	}
//#endif // PP_MATRIX_OUTPUT
//	cout << "Obj Function:\t";
//	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
//		cout << setw(PP_SETW) << PD_c[j];
//	if (PP_OUTPUT_LIMIT < PD_n)
//		cout << " ...";
//	cout << endl;
//	cout << "x0 =\t\t";
//	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_x0[j];
//	if (PP_OUTPUT_LIMIT < PD_n)
//		cout << " ...";
//	cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_x0);
//	cout << endl;
//	cout << "-------------------------------------------" << endl;
}

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PD_n; j++)
		parameterOutP->x[j] = parameterIn.x[j];
	parameterOutP->max_residual = parameterIn.max_residual;
}

// 0. Pseudoprojection
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx =";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter.x);

#ifdef PP_DEBUG
	cout << "\tActive:";
	for (int i = 0; i < PD_m; i++) {
		if (reduceResult->activeTag[i] == true)
			cout << " " << i;
	}
	cout << endl;
#else
	cout << endl;
#endif // PP_DEBUG
}

// 1. Movement on Polytope
void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
};

// 2. Least projection
void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx =";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter.x);

#ifdef PP_DEBUG
	cout << "\t\tActive:";
	for (int i = 0; i < PD_m; i++) {
		if (reduceResult->activeTag[i] == true)
			cout << " " << i;
	}
	cout << "\t";

	for (int i = 0; i < PD_m; i++)
		if (!PD_recessive_tag && PointInHalfspace(parameter.x, PD_A[i], PD_b[i])) {
			if (Vector_OnHyperplane(parameter.x, PD_A[i], PD_b[i]))
				continue;
			cout << "x inside half-space " << i << "!" << endl;
		}

	cout << "x belongs to hyperplanes: ";
	for (int i = 0; i < PD_n; i++) {
		if (PD_recessive_tag[i] && Vector_OnHyperplane(parameter.x, PD_A[i], PD_b[i]))
			cout << i << " ";
	}
	cout << endl;

#else
	cout << endl;
#endif // PP_DEBUG
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Start
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	if (PD_index % 1000 == 0)
		cout << PD_index << " problems processed. Time: " << (clock() - PD_time) / CLOCKS_PER_SEC << endl;
	//	cout << "Before problem output!" << endl;
	ProblemOutput(t);
}

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

// 2. Least projection
void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//---------------------------------- Problem functions -------------------------
inline PT_float_T Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
	PT_float_T sum = 0;
	for (int j = 0; j < PD_n; j++)
		sum += x[j] * y[j];
	return sum;
}

inline PT_float_T Vector_Norm(PT_vector_T x) {
	return sqrt(Vector_NormSquare(x));
}

inline PT_float_T Vector_NormSquare(PT_vector_T x) {
	PT_float_T sum = 0;

	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline bool PointInHalfspace // If the point belongs to the Halfspace with prescigion of PP_EPS_ZERO
(PT_vector_T x, PT_vector_T a, PT_float_T b) {
	return Vector_DotProduct(a, x) <= b + PP_EPS_ZERO;
}

inline bool PointInPolytope(PT_vector_T x) { // If the point belongs to the polytope with prescigion of PP_EPS_ZERO
	for (int i = 0; i < PD_m; i++) {
		if (PD_A[i][PP_ADD_FLAG] == 1)
			continue;
		if (!PointInHalfspace(x, PD_A[i], PD_b[i]))
			return false;
	}
	return true;
}

inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++)
		toPoint[j] = fromPoint[j];
}

inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PD_n; j++)
		equalVector[j] += plusVector[j];
}

inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PD_n; j++)
		equalPoint[j] -= minusVector[j];
}

inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] + y[j];
}

inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] - y[j];
}

inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] * r;
}

inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
	for (int j = 0; j < PD_n; j++)
		x[j] *= r;
}

inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
	for (int j = 0; j < PD_n; j++)
		x[j] /= r;
}

inline void Vector_ResetToZero(PT_vector_T x) {  // x = 0
	for (int j = 0; j < PD_n; j++) x[j] = 0;
}

inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] / r;
}

inline void Vector_Round(PT_vector_T x) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PD_n; j++) {
		if (fabs(x[j]) < PP_EPS_ZERO) {
			x[j] = 0;
			continue;
		}
		absValue = fabs(x[j]);
		sign = x[j] > 0 ? 1 : -1;
		floorValue = floor(absValue);
		fractionalPart = absValue - floorValue;
		if (1 - fractionalPart < PP_EPS_ZERO) {
			x[j] = sign * (floorValue + 1);
			continue;
		}
		if (fractionalPart < PP_EPS_ZERO)
			x[j] = sign * floorValue;
	}
}

inline void Vector_EpsZero(PT_vector_T x) { // If x[j] < PP_EPS_ZERO then x[j] = 0
	for (int j = 0; j < PD_n; j++)
		if (fabs(x[j]) < PP_EPS_ZERO)
			x[j] = 0;
}

inline PT_float_T ObjF(PT_vector_T x) {
	PT_float_T s = 0;
	for (int j = 0; j < PD_n; j++)
		s += PD_c[j] * x[j];
	return s;
}

//static bool LoadMatrixFormat() {
//	int nor,	// Number of matrix rows
//		noc,	// Number of matrix columns
//		non,	// Number of non-zero elements
//		noe;	// Number of equations
//	const char* mtxFile;
//	FILE* stream;// Input stream
//	char str[80] = { '\0' };
//	char* chr = str;
//
//	//--------------- Reading A ------------------
//	PD_MTX_File_A = PP_PATH;
//	PD_MTX_File_A += PP_MTX_PREFIX;
//	PD_MTX_File_A += PD_problemName;
//	PD_MTX_File_A += PP_MTX_POSTFIX_A;
//	mtxFile = PD_MTX_File_A.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Failure of opening file '" << mtxFile << "'.\n";
//		return false;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file " << mtxFile << endl;
//		return false;
//	}
//
//	if (nor >= noc) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "\n";
//		return false;
//	}
//
//	if (noc != PP_N) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Invalid input data: PP_N must be = " << noc << "\n";
//		return false;
//	}
//
//	if (nor != PP_M) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Invalid input data: PP_M must be = " << nor << "\n";
//		return false;
//	}
//
//	PD_m = noe = nor;
//	PD_n = noc;
//
//	if (2 * nor + noc > PP_MM) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Invalid input data: number of inequalities m = " << 2 * nor + noc
//			<< " must be < PP_MM + 1 =" << PP_MM + 1 << "\n";
//		return false;
//	}
//
//	for (int k = 0; k < non; k++) {
//		int i, j;
//
//		if (fscanf(stream, "%d%d%s", &i, &j, str) < 3) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file'" << mtxFile << "'." << endl;
//			return false;
//		}
//
//		i -= 1;
//		j -= 1;
//		if (i < 0) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Negative row index in'" << mtxFile << "'.\n" << endl;
//			return false;
//		}
//		if (j < 0) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Negative column index in'" << mtxFile << "'.\n" << endl;
//			return false;
//		}
//		PD_A[i][j] = strtod(str, &chr);
//	}
//
//	/*debug*
//for (int i = 0; i < PD_m; i++) {
//	for (int j = 0; j < PD_n; j++)
//		cout << PD_A[i][j] << " ";
//	cout << endl;
//}
///*end debug*/
//
//	fclose(stream);
//
//	//--------------- Reading b ------------------
//	PD_MTX_File_b = PP_PATH;
//	PD_MTX_File_b += PP_MTX_PREFIX;
//	PD_MTX_File_b += PD_problemName;
//	PD_MTX_File_b += PP_MTX_POSTFIX_B;
//	mtxFile = PD_MTX_File_b.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Failure of opening file '" << mtxFile << "'.\n";
//		return false;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file'" << mtxFile << "'." << endl;
//		return false;
//	}
//	if (noe != nor) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
//		return false;
//	}
//	if (noc != 1) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
//		return false;
//	}
//
//	for (int i = 0; i < noe; i++) {
//		if (fscanf(stream, "%s", str) < 1) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file '" << mtxFile << "'." << endl;
//			return false;
//		}
//		PD_b[i] = strtod(str, &chr);
//	}
//	fclose(stream);
//
//	/*debug*
//for (int i = 0; i < PD_m; i++)
//	cout << PD_b[i] << endl;
///*end debug*/
//
////--------------- Reading lo ------------------
//	PD_MTX_File_lo = PP_PATH;
//	PD_MTX_File_lo += PP_MTX_PREFIX;
//	PD_MTX_File_lo += PD_problemName;
//	PD_MTX_File_lo += PP_MTX_POSTFIX_LO;
//	mtxFile = PD_MTX_File_lo.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Failure of opening file '" << mtxFile << "'.\n";
//		return false;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file'" << mtxFile << "'." << endl;
//		return false;
//	}
//	if (nor != PD_n) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
//		return false;
//	}
//	if (noc != 1) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
//		return false;
//	}
//
//	for (int j = 0; j < PD_n; j++) {
//		if (fscanf(stream, "%s", str) < 1) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file '" << mtxFile << "'." << endl;
//			return false;
//		}
//		PD_lo[j] = strtod(str, &chr);
//	}
//
//	fclose(stream);
//
//	//--------------- Reading c ------------------
//	PD_MTX_File_c = PP_PATH;
//	PD_MTX_File_c += PP_MTX_PREFIX;
//	PD_MTX_File_c += PD_problemName;
//	PD_MTX_File_c += PP_MTX_POSTFIX_C;
//	mtxFile = PD_MTX_File_c.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Failure of opening file '" << mtxFile << "'.\n";
//		return false;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file'" << mtxFile << "'." << endl;
//		return false;
//	}
//	if (nor != PD_n) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
//		return false;
//	}
//	if (noc != 1) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
//		return false;
//	}
//
//	for (int j = 0; j < PD_n; j++) {
//		if (fscanf(stream, "%s", str) < 0) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file" << endl;
//			return false;
//		}
//		PD_c[j] = -strtod(str, &chr);
//	}
//	fclose(stream);
//
//	//--------------- Reading hi ------------------
//	PD_MTX_File_hi = PP_PATH;
//	PD_MTX_File_hi += PP_MTX_PREFIX;
//	PD_MTX_File_hi += PD_problemName;
//	PD_MTX_File_hi += PP_MTX_POSTFIX_HI;
//	mtxFile = PD_MTX_File_hi.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Failure of opening file '" << mtxFile << "'.\n";
//		return false;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file'" << mtxFile << "'." << endl;
//		return false;
//	}
//	if (nor != PD_n) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
//		return false;
//	}
//	if (noc != 1) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
//		return false;
//	}
//
//	for (int j = 0; j < PD_n; j++) {
//		if (fscanf(stream, "%s", str) < 1) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file '" << mtxFile << "'." << endl;
//			return false;
//		}
//		PD_hi[j] = strtod(str, &chr);
//	}
//	fclose(stream);
//
//	bool error = !Conversion();
//	if (error) return false;
//
//	//--------------- Reading x0 ------------------
//	PD_MTX_File_x0 = PP_PATH;
//	PD_MTX_File_x0 += PP_MTX_PREFIX;
//	PD_MTX_File_x0 += PD_problemName;
//	PD_MTX_File_x0 += PP_MTX_POSTFIX_X0;
//	mtxFile = PD_MTX_File_x0.c_str();
//	stream = fopen(mtxFile, "r+b");
//
//	if (stream == NULL) {
//		// Generating Coordinates of starting point
//		for (int j = 0; j < PD_n; j++)
//			PD_x0[j] = 0;
//		return true;
//	}
//
//	SkipComments(stream);
//	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Unexpected end of file'" << mtxFile << "'." << endl;
//		return false;
//	}
//	if (nor != PD_n) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
//		return false;
//	}
//	if (noc != 1) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout
//			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
//		return false;
//	}
//
//	for (int j = 0; j < PD_n; j++) {
//		if (fscanf(stream, "%s", str) < 0) {
//			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//				cout
//				<< "Unexpected end of file" << endl;
//			return false;
//		}
//		PD_x0[j] = strtod(str, &chr);
//	}
//	fclose(stream);
//	return true;
//}

static bool Conversion() { // Transformation to inequalities & dimensionality reduction
	static PT_float_T fvA[PP_MM]; // Free variable coefficients
	static bool Flag[PP_N];		// Flags of free variables to delete
	static int fvEqI;	// Inequality index of free variable
	static bool single;

	for (int jc = 0; jc < PD_n; jc++)
		if (PD_c[jc] == 0 && PD_hi[jc] == PP_INFINITY && PD_lo[jc] == 0)
			Flag[jc] = true;

	for (int jc = 0; jc < PD_n; jc++) {
		if (!Flag[jc])
			continue;
		for (int i = 0; i < PD_m; i++) { // Find corresponding equation
			if (PD_A[i][jc] == 0)
				continue;
			single = true;
			for (int ii = i + 1; ii < PD_m; ii++) // Vertical uniqueness
				if (PD_A[ii][jc] != 0) {
					single = false;
					break;
				}
			if (single)
				fvEqI = i;
			break;
		}

		if (!single) {
			Flag[jc] = false;
		}
		else
			fvA[fvEqI] = PD_A[fvEqI][jc];
	}

	static bool PD_delete[PP_MM]; // Rows to delete
	PT_float_T s;

	for (int i = 0; i < PD_m; i++) { // Check inconsistent end degenerate equation
		s = 0;
		for (int j = 0; j < PD_n; j++)
			s += fabs(PD_A[i][j]);
		if (s == 0) {
			if (PD_b[i] != 0) {
				//
				cout
					<< "Inconsistent equation " << i << ": " << s << " = " << PD_b[i] << endl;
				return false;
				PD_delete[i] = true;
			}
		}
	}

	for (int i = 0; i < PD_m; i++) { // Removing degenerate equations
		if (!PD_delete[i]) continue;
		for (int ii = i; ii < PD_m - 1; ii++) {  // Remove null equation
			for (int j = 0; j < PD_n; j++)
				PD_A[ii][j] = PD_A[ii + 1][j];
			fvA[ii] = fvA[ii + 1];
			PD_b[ii] = PD_b[ii + 1];
		}
		PD_m--;
	}

	for (int jc = 0; jc < PD_n; jc++) { // Delete free variables
		if (!Flag[jc])
			continue;
		for (int j = jc; j < PD_n - 1; j++) { // Delete column
			PD_c[j] = PD_c[j + 1];
			PD_lo[j] = PD_lo[j + 1];
			PD_hi[j] = PD_hi[j + 1];
			Flag[j] = Flag[j + 1];
			for (int i = 0; i < PD_m; i++)
				PD_A[i][j] = PD_A[i][j + 1];
		}

		PD_n--;
		jc--;
		for (int i = 0; i < PD_m; i++)
			PD_A[i][PD_n] = 0;
		PD_c[PD_n] = 0;
		PD_lo[PD_n] = 0;
		PD_hi[PD_n] = 0;
	}

	int m = PD_m;
	for (int i = 0; i < m; i++) { // Conversion to inequalities

		if (fvA[i] == 0) { // Equation without free variable => adding inequality.
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = -PD_A[i][j];
			PD_b[PD_m] = -PD_b[i];
			PD_m++;
			assert(PD_m <= PP_MM);
		}
		else {
			if (fvA[i] < 0) { // Free variable is negative => change sign to opposite.
				for (int j = 0; j < PD_n; j++)
					PD_A[i][j] = -PD_A[i][j];
				PD_b[i] = -PD_b[i];
			}
		}
	}

	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
		for (int j = 0; j < PD_n; j++)
			if (PD_A[i][j] == 0)
				PD_A[i][j] = 0;

	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
		for (int j = 0; j < PD_n; j++)
			PD_A[i + PD_m][j] = 0;
		PD_A[i + PD_m][i] = -1;
		if (PD_lo[i] == 0)
			PD_b[i + PD_m] = 0;
		else
			PD_b[i + PD_m] = -PD_lo[i];
	}
	PD_m += PD_n;
	assert(PD_m <= PP_MM);

	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
		if (PD_hi[i] != PP_INFINITY) {
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = 0;
			PD_A[PD_m][i] = 1;
			PD_b[PD_m] = PD_hi[i];
			PD_m++;
			assert(PD_m <= PP_MM);
		}
	}
	return true;
}

//inline bool MTX_SavePoint(PT_vector_T x, const char* filename, const char* comment) {
//	FILE* stream;
//
//	stream = fopen(filename, "w");
//	if (stream == NULL) {
//		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
//			cout << "Failure of opening file '" << filename << "'.\n";
//		return false;
//	}
//
//	fprintf(stream, "%c %s\n%d %d\n", '%', comment, PD_n, 1);
//
//	for (int j = 0; j < PD_n; j++)
//		fprintf(stream, "%.16f\n", x[j]);
//
//	fclose(stream);
//	return true;
//}

inline bool Vector_OnHyperplane // If the point belongs to the Hyperplane with prescigion of PP_EPS_ZERO
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return fabs(Vector_DotProduct(a, point) - b) < PP_EPS_ZERO;
}

// Vector r of orthogonal projection of point z onto Half-space <a,x> <= b
inline PT_float_T // residual
Vector_OrthogonalProjectionOntoHalfspace(PT_vector_T z, PT_vector_T a, PT_float_T b, PT_vector_T r, int* exitCode) {
	PT_float_T factor;
	PT_float_T aNormSquare = Vector_NormSquare(a); // ||a||^2
	PT_float_T residual = Vector_DotProduct(a, z) - b; // <a,z> - b

	if (sqrt(aNormSquare) < PP_EPS_ZERO) {
		*exitCode = PP_EXITCODE_DEGENERATE_INEQUALITY;
		for (int j = 0; j < PD_n; j++)
			r[j] = 0;
		return -1;
	}

	if (residual < 0) { // <a,z> - b < 0
		*exitCode = PP_EXITCODE_INSIDE_HALFSPACE;
		for (int j = 0; j < PD_n; j++)
			r[j] = 0;
		return residual;
	}

	factor = -residual / aNormSquare; // (b - <z,a>) / ||a||^2
	for (int j = 0; j < PD_n; j++) // a(b - <z,a>) / ||a||^2
		r[j] = factor * a[j];
	*exitCode = PP_EXITCODE_SUCCESSFUL_PROJECTING;
	return residual;
}

// Distance from point z to halfspace <a,x> <= b: |<a,z> - b|/||a||
inline PT_float_T Vector_DistanceToHalfspace(PT_vector_T z, PT_vector_T a, PT_float_T b) {
	PT_float_T aNorm = sqrt(Vector_NormSquare(a)); // ||a||
	PT_float_T a_dot_z_minus_b; // <a,z> - b

	if (aNorm < PP_EPS_ZERO) //Degenerate equation
		return 0;

	a_dot_z_minus_b = Vector_DotProduct(a, z) - b; // <a,z> - b
	if (a_dot_z_minus_b <= PP_EPS_ZERO) // Point belongs to halfspace
		return 0;

	return a_dot_z_minus_b / aNorm;
}

// Vector o of oblique projection of point z onto Half-space <a,x> <= b with respect to vector g
inline void Vector_ObliqueProjectionOntoHalfspace(PT_vector_T z, PT_vector_T a, PT_float_T b, PT_vector_T g, PT_vector_T o, int* exitCode) {
	PT_float_T a_dot_g;	// <a,g>
	PT_float_T a_dot_z_minus_b;	// <a,z> - b
	PT_float_T factor;	// (b - <a,z>) / <a,g>

	a_dot_z_minus_b = Vector_DotProduct(a, z) - b; // <a,z> - b
	if (a_dot_z_minus_b <= -PP_EPS_ZERO) { // <a,z> - b < 0
		*exitCode = PP_EXITCODE_INSIDE_HALFSPACE;
		for (int j = 0; j < PD_n; j++)
			o[j] = 0;
		return;
	}

	if (fabs(a_dot_z_minus_b) < PP_EPS_ZERO) { // <a,z> - b < 0
		*exitCode = PP_EXITCODE_ON_HYPERPLANE;
		for (int j = 0; j < PD_n; j++)
			o[j] = 0;
		return;
	}

	a_dot_g = Vector_DotProduct(a, g); // <a,g>
	if (fabs(a_dot_g) < PP_EPS_ZERO) {
		*exitCode = PP_EXITCODE_PARALLEL;
		for (int j = 0; j < PD_n; j++)
			o[j] = PP_INFINITY;
		return;
	}

	if (a_dot_g >= PP_EPS_ZERO) {
		*exitCode = PP_EXITCODE_RECESSIVE;
		for (int j = 0; j < PD_n; j++)
			o[j] = PP_INFINITY;
		return;
	}

	factor = a_dot_z_minus_b / a_dot_g; // (<a,z> - b) / <a,g>

	// Oblique projection vector: o = -(<a,z> - b)g/<a, g> = -factor * r
	for (int j = 0; j < PD_n; j++)
		o[j] = -factor * g[j];

	*exitCode = PP_EXITCODE_SUCCESSFUL_PROJECTING;
	return;
};

inline PT_float_T Distance(PT_vector_T x, PT_vector_T y) {
	PT_vector_T z;
	Vector_Subtraction(x, y, z);
	return Vector_Norm(z);
}

inline void UnitObjVector(PT_vector_T objUnitVector) { // Calculating Objective Unit Vector
	double c_norm = Vector_Norm(PD_c);
	Vector_DivideByNumber(PD_c, c_norm, objUnitVector);
}

inline void ProblemOutput(double elapsedTime) {
	Vector_Round(PD_u);
	
	//cout << "=============================================" << endl;
	//cout << "Elapsed time: " << elapsedTime << endl;
	//cout << "Iterations: " << BSF_sv_iterCounter << endl;
	//cout << "=============================================" << endl;

	//PD_MTX_File_sp = PP_PATH;
	//PD_MTX_File_sp += PP_MTX_PREFIX;
	//PD_MTX_File_sp += PD_problemName;
	//PD_MTX_File_sp += PP_MTX_POSTFIX_U0;
	//const char* sufacePointFile = PD_MTX_File_sp.c_str();
	//const char* comment = "Surface point";
	//if (MTX_SavePoint(PD_u, sufacePointFile, comment))
	//	cout << "Surface point is saved into the file '" << sufacePointFile << "'." << endl;
#ifdef PP_DATABASE_OUTPUT
	SurfacePoint point;
	vector<double> _buff(PD_n);

	for (unsigned i = 0; i < PD_n; i++)
		_buff[i] = PD_u[i];
	point.coefficients = doubleToChar(_buff);
	point.problem_id = PD_ids[PD_index - 1];

	point.id = storage.insert(point);
	//if (point.id == 0)
	//	cout << "Error: new surface point has not been inserted into SQLite DB!" << endl;
#else
	PD_currentProblem->x0 = new CArray(PD_n);
	for (unsigned i = 0; i < PD_n; i++)
		PD_currentProblem->x0->setValue(i, PD_u[i]);
	PD_packetWriter->addProblem(*PD_currentProblem);
	PD_currentProblem->x0->~CArray();
#endif // PP_DATABASE_OUTPUT


	//cout << "Surface point u0:\t";
	//for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_u[j];
	//if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	//cout << endl;
	//cout << "Polytope residual: " << PolytopeResidual(PD_u) << endl;

#ifndef PP_DATABASE_OUTPUT
	if (PD_index == PD_packetReader->packetSize)
		PD_packetWriter->close();
#endif // !PP_DATABASE_OUTPUT

}

//inline void SkipComments(FILE* stream) {
//	fpos_t pos;	// Position in the input stream
//	int res;
//	res = fscanf(stream, "\n");
//	fgetpos(stream, &pos);
//	while (getc(stream) == '%') {
//		while (getc(stream) != 10);
//		res = fscanf(stream, "\n");
//		fgetpos(stream, &pos);
//	};
//	fsetpos(stream, &pos);
//}

inline PT_float_T PolytopeResidual(PT_vector_T x) { // Measure of distance from point to polytope
	PT_float_T sum = 0;
	int nonzero = 0;
	PT_float_T distance;

	for (int i = 0; i < PD_m; i++) {
		distance = Vector_DistanceToHalfspace(x, PD_A[i], PD_b[i]);
		if (distance > 0) {
			sum += distance;
			nonzero++;
		}
	}
	if (nonzero != 0)
		sum /= nonzero;
	return sum;
}

inline void ApexPoint(PT_vector_T innerPont, PT_vector_T apexPoint) {
	PT_float_T a_dot_e_c, a_dot_innerPoint;
	PT_float_T max_cDistance = 0;
	PT_float_T cFactor;

	for (int i = 0; i < PD_m; i++) {
		a_dot_e_c = Vector_DotProduct(PD_A[i], PD_e_c);
		if (a_dot_e_c <= -PP_EPS_ZERO) {// not recessive!
			PD_recessive_tag[i] = false;
			continue;
		}
		PD_recessive_tag[i] = true;
		a_dot_innerPoint = Vector_DotProduct(PD_A[i], innerPont);
		cFactor = (PD_b[i] - a_dot_innerPoint) / a_dot_e_c;
		assert(cFactor > -(PP_EPS_ZERO * 10));
		max_cDistance = PF_MAX(max_cDistance, cFactor);
	}
	max_cDistance += PP_ETA_TO_APEX;
	Vector_MultiplyByNumber(PD_e_c, max_cDistance, PD_direction);
	Vector_Addition(innerPont, PD_direction, apexPoint);
}

#ifdef PP_DATABASE_OUTPUT
std::vector<double> charToDouble(std::vector<char> _In) {
	std::vector<double> _Out;
	size_t size = _In.size();
	_Out.resize(size / sizeof(double));
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

std::vector<char> doubleToChar(std::vector<double> _In) {
	std::vector<char> _Out;
	size_t size = _In.size() * sizeof(double);
	_Out.resize(size);
	std::memcpy(_Out.data(), _In.data(), size);
	return _Out;
}

void printLppForm(Problem problem, std::vector<Inequality> inequalities) {
	int M = 2 * problem.N + inequalities.size();
	int width = 10;
	std::vector<double> _vec;
	std::cout << problem.id << '\t' << problem.N << '\t' << M << std::endl;
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.high << std::endl;
	}
	for (int i = 0; i < inequalities.size(); i++) {
		_vec = charToDouble(inequalities[i].coefficients);
		for (int j = 0; j < problem.N; j++)
			std::cout << std::setw(width) << _vec[j];
		std::cout << std::setw(width) << inequalities[i].b << std::endl;
	}
	for (int i = 0; i < problem.N; i++) {
		for (int j = 0; j < problem.N; j++)
			if (i == j) std::cout << std::setw(width) << double(-1);
			else std::cout << std::setw(width) << double(0);
		std::cout << std::setw(width) << problem.low << std::endl;
	}
	_vec = charToDouble(problem.c);
	std::copy(_vec.begin(), _vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
	std::cout << std::endl;
}
#endif // PP_DATABASE_OUTPUT
// MPI
// MPI
// MPI
// MPI
// MPI
// MPI