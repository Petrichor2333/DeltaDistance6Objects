#include "pch.h"
#include "config.h"
#include "SolvingMinimumMove.h"
#include "SolvingMinimum.h"
#include "DDFunc.h"

void ValueOfVCCD(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* val, double R[][3], double* p0, double* v, double t) {
	double p[3];
	for (int i = 0;i < 3;i++) {
		p[i] = p0[i] + t * v[i];
	}

	double x_in_2[3] = { 0,0,0 };
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			x_in_2[i] += R[i][j] * x[j];
		}
		x_in_2[i] += p[i];
	}
	double Val1;
	obj1->DDFuncSqValu(x, &Val1);
//	DDFuncSqValu_Obj1_MS(x, &Val1);
//	DDFuncSqValu_Obj2_PC(x, &Val1);
//	DDFuncSqValu_Obj3_DM(x, &Val1);
//	DDFuncSqValu_Obj4_CK(x, &Val1);
//	DDFuncSqValu_Obj5_PL(x, &Val1);
//	DDFuncSqValu_Obj6_HN(x, &Val1);

	double Val2;
	obj2->DDFuncSqValu(x_in_2, &Val2);
//	DDFuncSqValu_Obj1_MS(x_in_2, &Val2);
//	DDFuncSqValu_Obj2_PC(x_in_2, &Val2);
//	DDFuncSqValu_Obj3_DM(x_in_2, &Val2);
//	DDFuncSqValu_Obj4_CK(x_in_2, &Val2);
//	DDFuncSqValu_Obj5_PL(x_in_2, &Val2);
//	DDFuncSqValu_Obj6_HN(x_in_2, &Val2);
	*val = Val1 + Val2;
	return;
}

void GradOfVCCD(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* HGrad, double R[][3], double* p0, double* v, double t) {
	double p[3];
	for (int i = 0;i < 3;i++) {
		p[i] = p0[i] + t * v[i];
	}

	double x_in_2[3] = { 0,0,0 };
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			x_in_2[i] += R[i][j] * x[j];
		}
		x_in_2[i] += p[i];
	}

	double Grad1[3], Grad2[3];

	obj1->DDFuncSqGrad(x, Grad1);
//	DDFuncSqGrad_Obj1_MS(x, Grad1);
//	DDFuncSqGrad_Obj2_PC(x, Grad1);
//	DDFuncSqGrad_Obj3_DM(x, Grad1);
//	DDFuncSqGrad_Obj4_CK(x, Grad1);
//	DDFuncSqGrad_Obj5_PL(x, Grad1);
//	DDFuncSqGrad_Obj6_HN(x, Grad1);

	obj2->DDFuncSqGrad(x_in_2, Grad2);
//	DDFuncSqGrad_Obj1_MS(x_in_2, Grad2);
//	DDFuncSqGrad_Obj2_PC(x_in_2, Grad2);
//	DDFuncSqGrad_Obj3_DM(x_in_2, Grad2);
//	DDFuncSqGrad_Obj4_CK(x_in_2, Grad2);
//	DDFuncSqGrad_Obj5_PL(x_in_2, Grad2);
//	DDFuncSqGrad_Obj6_HN(x_in_2, Grad2);

	for (int i = 0;i < 3;i++) {
		HGrad[i] = Grad1[i];
		for (int j = 0;j < 3;j++) {
			HGrad[i] += R[j][i] * Grad2[j];
		}
	}
	HGrad[3] = 0;
	for (int i = 0;i < 3;i++) {
		HGrad[3] += Grad2[i] * v[i];
	}

	return;
}

void HessOfVCCD(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* Hess, double R[][3], double* p0, double* v, double t) {
	double p[3];
	for (int i = 0;i < 3;i++) {
		p[i] = p0[i] + t * v[i];
	}
	
	double x_in_2[3] = { 0,0,0 };
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			x_in_2[i] += R[i][j] * x[j];
		}
		x_in_2[i] += p[i];
	}

	double Hess1[6], Hess2[6];

	obj1->DDFuncSqHess(x, Hess1);
//	DDFuncSqHess_Obj1_MS(x, Hess1);
//	DDFuncSqHess_Obj2_PC(x, Hess1);
//	DDFuncSqHess_Obj3_DM(x, Hess1);
//	DDFuncSqHess_Obj4_CK(x, Hess1);
//	DDFuncSqHess_Obj5_PL(x, Hess1);
//	DDFuncSqHess_Obj6_HN(x, Hess1);

	obj2->DDFuncSqHess(x_in_2, Hess2);
//	DDFuncSqHess_Obj1_MS(x_in_2, Hess2);
//	DDFuncSqHess_Obj2_PC(x_in_2, Hess2);
//	DDFuncSqHess_Obj3_DM(x_in_2, Hess2);
//	DDFuncSqHess_Obj4_CK(x_in_2, Hess2);
//	DDFuncSqHess_Obj5_PL(x_in_2, Hess2);
//	DDFuncSqHess_Obj6_HN(x_in_2, Hess2);

	double RTH2[9];
	RTH2[0] = Hess2[0] * R[0][0] + Hess2[1] * R[1][0] + Hess2[2] * R[2][0];
	RTH2[1] = Hess2[1] * R[0][0] + Hess2[3] * R[1][0] + Hess2[4] * R[2][0];
	RTH2[2] = Hess2[2] * R[0][0] + Hess2[4] * R[1][0] + Hess2[5] * R[2][0];
	RTH2[3] = Hess2[0] * R[0][1] + Hess2[1] * R[1][1] + Hess2[2] * R[2][1];
	RTH2[4] = Hess2[1] * R[0][1] + Hess2[3] * R[1][1] + Hess2[4] * R[2][1];
	RTH2[5] = Hess2[2] * R[0][1] + Hess2[4] * R[1][1] + Hess2[5] * R[2][1];
	RTH2[6] = Hess2[0] * R[0][2] + Hess2[1] * R[1][2] + Hess2[2] * R[2][2];
	RTH2[7] = Hess2[1] * R[0][2] + Hess2[3] * R[1][2] + Hess2[4] * R[2][2];
	RTH2[8] = Hess2[2] * R[0][2] + Hess2[4] * R[1][2] + Hess2[5] * R[2][2];

	Hess[0] = Hess1[0] + R[0][0] * RTH2[0] + R[1][0] * RTH2[1] + R[2][0] * RTH2[2];
	Hess[1] = Hess1[1] + R[0][1] * RTH2[0] + R[1][1] * RTH2[1] + R[2][1] * RTH2[2];
	Hess[2] = Hess1[2] + R[0][2] * RTH2[0] + R[1][2] * RTH2[1] + R[2][2] * RTH2[2];
	Hess[4] = Hess1[3] + R[0][1] * RTH2[3] + R[1][1] * RTH2[4] + R[2][1] * RTH2[5];
	Hess[5] = Hess1[4] + R[0][2] * RTH2[3] + R[1][2] * RTH2[4] + R[2][2] * RTH2[5];
	Hess[7] = Hess1[5] + R[0][2] * RTH2[6] + R[1][2] * RTH2[7] + R[2][2] * RTH2[8];

	double H2v[3];
	H2v[0] = Hess2[0] * v[0] + Hess2[1] * v[1] + Hess2[2] * v[2];
	H2v[1] = Hess2[1] * v[0] + Hess2[3] * v[1] + Hess2[4] * v[2];
	H2v[2] = Hess2[2] * v[0] + Hess2[4] * v[1] + Hess2[5] * v[2];

	Hess[3] = 0, Hess[6] = 0, Hess[8] = 0;
	for (int j = 0;j < 3;j++) {
		Hess[3] += R[j][0] * H2v[j];
		Hess[6] += R[j][1] * H2v[j];
		Hess[8] += R[j][2] * H2v[j];
	}

	Hess[9] = v[0] * H2v[0] + v[1] * H2v[1] + v[2] * H2v[2];

	return;
}

bool CholeskySolveNewtonStepCCD(double* H, double* g, double* x) {
	double L[10];
	double l0_sq, l4_sq, l7_sq, l9_sq;

	l0_sq = H[0];
	if (l0_sq < fabs(l0_sq)*epsilon) { return 0; }
	L[0] = sqrt(l0_sq);
	L[1] = H[1] / L[0];
	L[2] = H[2] / L[0];
	L[3] = H[3] / L[0];
	l4_sq = H[4] - L[1] * L[1];
	if (l4_sq < fabs(l4_sq)*epsilon) { return 0; }
	L[4] = sqrt(l4_sq);
	L[5] = (H[5] - L[1] * L[2]) / L[4];
	L[6] = (H[6] - L[1] * L[3]) / L[4];
	l7_sq = H[7] - L[2] * L[2] - L[5] * L[5];
	if (l7_sq < fabs(l7_sq)*epsilon) { return 0; }
	L[7] = sqrt(l7_sq);
	L[8] = (H[8] - L[2] * L[3] - L[5] * L[6]) / L[7];
	l9_sq = H[9] - L[3] * L[3] - L[6] * L[6] - L[8] * L[8];
	if (l9_sq < fabs(l9_sq)*epsilon) { return 0; }
	L[9] = sqrt(l9_sq);

	double y[4];
	y[0] = -g[0] / L[0];										//x=-Hess^{-1}Grad, so "-g"
	y[1] = (-g[1] - L[1] * y[0]) / L[4];
	y[2] = (-g[2] - L[2] * y[0] - L[5] * y[1]) / L[7];
	y[3] = (-g[3] - L[3] * y[0] - L[6] * y[1] - L[8] * y[2]) / L[9];

	x[3] = y[3] / L[9];
	x[2] = (y[2] - L[8] * x[3]) / L[7];
	x[1] = (y[1] - L[5] * x[2] - L[6] * x[3]) / L[4];
	x[0] = (y[0] - L[1] * x[1] - L[2] * x[2] - L[3] * x[3]) / L[0];
	return 1;
}

bool SolveNewtonStepCCD(double* H, double* g, double* x) {

	if (CholeskySolveNewtonStepCCD(H, g, x)) {
		return 1;
	}

	integer mat_size = 4;
	char JOBZ = 'V';
	char UPLO = 'U';
	double P[16] = {H[0], H[1], H[2], H[3],
					H[1], H[4], H[5], H[6],
					H[2], H[5], H[7], H[8] ,
					H[3], H[6], H[8], H[9]};
	integer LDA = 4;
	double eig_val[4];
	double WORK[16];
	integer LWORK = 16;
	integer Info;
	int info = 0;

	info = dsyev_(&JOBZ, &UPLO, &mat_size, P, &LDA, eig_val, WORK, &LWORK, &Info);
	if (info != 0) {

		cout << "Eigen Decomposition Solve Newton Step Failed..........." << endl;
		cout << endl;
		//assert(0);

		return 0;
	}

	//row of P is eigen vector, H=P^{T}*D*P
	double y[4], z[4];
	for (int i = 0;i < 4;i++) {
		y[i] = 0;
		for (int j = 0;j < 4;j++) {
			y[i] -= P[4 * i + j] * g[j];
		}
	}

	//	double tolerance = epsilon * fabs(eig_val[2]);
	for (int i = 0;i < 4;i++) {
		z[i] = eig_val[i] > epsilon * fabs(eig_val[i]) ? y[i] / fabs(eig_val[i]) : 0;
	}

	for (int i = 0;i < 4;i++) {
		x[i] = 0;
		for (int j = 0;j < 4;j++) {
			x[i] += P[4 * j + i] * z[j];
		}
	}
	return 1;
}

bool SolveForCollisionMove(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* xInit, ResultsOfSolvingMinimumMove * Result, double R1[][3], double* p10, double R2[][3], double* p20, double*v1, double* v2) {

	double R[3][3];		//frame 1 in 2
	double p0[3];
	double v[3];

	for (int i = 0;i < 3;i++) {
		p0[i] = 0;
		v[i] = 0;
		for (int j = 0;j < 3;j++) {
			p0[i] += R2[j][i] * (p10[j] - p20[j]);
			v[i] += R2[j][i] * (v1[j] - v2[j]);
			R[i][j] = 0;
			for (int k = 0;k < 3;k++) {
				R[i][j] += R2[k][i] * R1[k][j];
			}
		}

	}

	double x[3];
	for (int i = 0;i < 3;i++) {
		x[i] = 0;
		for (int j = 0;j < 3;j++) {
			x[i] += R1[j][i] * (xInit[j] - p10[j]);
		}
	}

	double t = 0;

	double HGrad[4];	//Homogeneous Gradient
	GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);
	double GradNorm = sqrt(HGrad[0] * HGrad[0] + HGrad[1] * HGrad[1] + HGrad[2] * HGrad[2] + HGrad[3] * HGrad[3]);

	double q[4];	//Newton Step

	int icount = 0;

	while (GradNorm > ERROR_GRAD && icount < 101) {
		double HHess[10];	//Homogeneous Hessian matrix
		HessOfVCCD(obj1, obj2, x, HHess, R, p0, v, t);
		if (!SolveNewtonStepCCD(HHess, HGrad, q)) {		//Solve the Newton step p=-Hess^{+}Grad
			Result->if_success = 0;
			return 0;
		}

		double ModPara = 0.1;

		double q_norm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		double max_step = 5;
		if (q_norm > max_step) {
			q[0] *= (max_step / q_norm);
			q[1] *= (max_step / q_norm);
			q[2] *= (max_step / q_norm);
			q[3] *= (max_step / q_norm);
			q_norm = max_step;
		}

		double V_x;
		ValueOfVCCD(obj1, obj2, x, &V_x, R, p0, v, t);
		double alpha = 1, ro = 0.8;
		double V_al;
		double x_al[3], t_al;
		x_al[0] = x[0] + alpha * q[0], x_al[1] = x[1] + alpha * q[1], x_al[2] = x[2] + alpha * q[2], t_al = t + alpha * q[3];
		ValueOfVCCD(obj1, obj2, x_al, &V_al, R, p0, v, t_al);
		double  alpha_b = ro * alpha;
		double V_al_b;
		double x_al_b[3], t_al_b;
		x_al_b[0] = x[0] + alpha_b * q[0], x_al_b[1] = x[1] + alpha_b * q[1], x_al_b[2] = x[2] + alpha_b * q[2], t_al_b = t + alpha_b * q[3];
		ValueOfVCCD(obj1, obj2, x_al_b, &V_al_b, R, p0, v, t_al_b);
		while (V_al - V_al_b > 0 && alpha > 0.05) {
			alpha = alpha_b;
			x_al[0] = x_al_b[0], x_al[1] = x_al_b[1], x_al[2] = x_al_b[2], t_al = t_al_b;
			V_al = V_al_b;
			alpha_b *= ro;
			x_al_b[0] = x[0] + alpha_b * q[0], x_al_b[1] = x[1] + alpha_b * q[1], x_al_b[2] = x[2] + alpha_b * q[2], t_al_b = t + alpha_b * q[3];
			ValueOfVCCD(obj1, obj2, x_al_b, &V_al_b, R, p0, v, t_al_b);
		}
		x[0] = x_al[0], x[1] = x_al[1], x[2] = x_al[2], t = t_al;

		GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);
		GradNorm = sqrt(HGrad[0] * HGrad[0] + HGrad[1] * HGrad[1] + HGrad[2] * HGrad[2] + HGrad[3] * HGrad[3]);

		icount++;
	}
	if (icount > 49) {
		Result->if_success = 0;
		cout << "--------------------------------------------------" << endl;
		cout << "Failed: icount = " << icount << " > 49" << endl;
		cout << "Failed V: --------------" << endl;
		double V_fail;
		ValueOfVCCD(obj1, obj2, x, &V_fail, R, p0, v, t);
		cout << V_fail << endl;
		cout << endl;
		cout << "Failed HGrad:" << endl;
		cout << HGrad[0] << ", " << HGrad[1] << ", " << HGrad[2] << ", " << HGrad[3] << endl;
		cout << endl;
		cout << "Failed q:" << endl;
		cout << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << endl;
		cout << endl;
		cout << "x:" << endl;
		cout << x[0] << ", " << x[1] << ", " << x[2] << endl;
		cout << "t:" << endl;
		cout << t << endl;
		cout << "--------------------------------------------------" << endl;
	}
	else {
		Result->if_success = 1;
	}
	ValueOfVCCD(obj1, obj2, x, &(Result->V_min), R, p0, v, t);
	if (Result->V_min < ERROR_COLLIDE_V && t >= 0) {
		Result->if_collide = 1;
	}
	else {
		Result->if_collide = 0;
	}
	Result->num_iter = icount;
	for (int i = 0;i < 3;i++) {
		xInit[i] = p0[i] + t * v[i];
		for (int j = 0;j < 3;j++) {
			xInit[i] += R1[i][j] * x[j];
		}
	}

	if (Result->if_collide != 1) {
		return 1;
	}



	return 1;
}