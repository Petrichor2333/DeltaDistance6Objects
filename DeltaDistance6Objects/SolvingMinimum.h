#pragma once

#include "DDFunc.h"

void ValueOfV(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* val, double R[][3], double* t);
void GradOfV(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* Grad, double R[][3], double* t);
void HessOfV(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* Hess, double R[][3], double* t);

bool CholeskySolveNewtonStep(double* H, double* g, double* x);
bool SolveNewtonStep(double* H, double* g, double* x);

struct ResultsOfSolvingMinimum {
	bool if_success;
	double ConvPoint[3];
	double V_min;
	bool if_collide;
	int num_iter;
};

// Solving the problem of collision detection via unconstrained optimization
bool SolveForCollision(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* xInit, ResultsOfSolvingMinimum * Result, double R1[][3], double* t1, double R2[][3], double* t2);



