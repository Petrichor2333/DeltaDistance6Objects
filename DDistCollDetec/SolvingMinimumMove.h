#pragma once

#include "DDFunc.h"

//void ValueOfVatXMove(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* val, double R[][3], double* p0, double* v, double t);

bool CholeskySolveNewtonStepMove(double* H, double* g, double* x);
bool EigenDecompSolveNewtonStepMove(double* H, double* g, double* x);

struct ResultsOfSolvingMinimumMove {
	bool if_success;
	double* ConvPoint;
	double V_min;
	bool if_collide;
	int num_iter;
};

// Solving the problem of collision detection for moving objects via unconstrained optimization
bool SolveForCollisionMove(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, ResultsOfSolvingMinimumMove * Result, double R1[][3], double* p10, double R2[][3], double* p2, double*v1, double* v2);



