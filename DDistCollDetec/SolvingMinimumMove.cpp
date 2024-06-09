#include "pch.h"
#include "config.h"
#include "SolvingMinimumMove.h"
#include "SolvingMinimum.h"
#include "DDFunc.h"

/*
void ValueOfVatXMove(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* val, double R[][3], double* p0, double* v, double t) {
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
	obj1->DDFuncSqValuIndep(x, &Val1, obj1->PObjData);

	double Val2;
	obj2->DDFuncSqValuIndep(x_in_2, &Val2, obj2->PObjData);

	*val = Val1 + Val2;
	return;
}
*/

bool CholeskySolveNewtonStepMove(double* H, double* g, double* x) {
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

bool EigenDecompSolveNewtonStepMove(double* H, double* g, double* x) {
	integer mat_size = 4;
	char JOBZ = 'V';
	char UPLO = 'U';
	double P[16] = { H[0], H[1], H[2], H[3],
					H[1], H[4], H[5], H[6],
					H[2], H[5], H[7], H[8] ,
					H[3], H[6], H[8], H[9] };
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

bool SolveForCollisionMove(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, ResultsOfSolvingMinimumMove * Result, double R1[][3], double* p10, double R2[][3], double* p20, double*v1, double* v2) {

	double R[3][3];		//frame 1 in 2
	double p0[3];
	double v[3];   //velocity of 1 wrt 2, represented in 2
	double v2_in_2[3];	//velocity of 2 wrt 1, represented in 2;
	double x[3];

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
		v2_in_2[i] = -v[i];
	}


	for (int i = 0;i < 3;i++) {
		x[i] = 0;
		for (int j = 0;j < 3;j++) {
			x[i] -= R[j][i] * p0[j];
		}
		x[i] *= 0.5;
	}
	double t = 0;

	/*
	double p[3];
	for (int i = 0;i < 3;i++) {
		p[i] = p0[i] + t * v[i];
	}*/
	
	double x_in_2[3] = { 0,0,0 };
	for (int i = 0;i < 3;i++) {
		for (int j = 0;j < 3;j++) {
			x_in_2[i] += R[i][j] * x[j];
		}
		x_in_2[i] += p0[i];
	}

	double Val1, Val2;

	double HGrad[4];	//Homogeneous Gradient
	double Grad1[3], Grad2[3];
	//GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);
	/*-------------------------GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);------------------------*/
	obj1->DDFuncSqGrad(x, Grad1, obj1->PObjData);
	obj2->DDFuncSqGrad(x_in_2, Grad2, obj2->PObjData);

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
	/*-------------------------GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);------------------------*/

	double GradNorm = sqrt(HGrad[0] * HGrad[0] + HGrad[1] * HGrad[1] + HGrad[2] * HGrad[2] + HGrad[3] * HGrad[3]);

	double q[4];	//Newton Step

	int icount = 0;

	while (GradNorm > ERROR_GRAD && icount < 55) {
		double HHess[10];	//Homogeneous Hessian matrix
		
		//HessOfVCCD(obj1, obj2, x, HHess, R, p0, v, t);
		/*-------------------------HessOfVCCD(obj1, obj2, x, HHess, R, p0, v, t);------------------------*/
		double Hess1[6], Hess2[6];
		obj1->DDFuncSqHess(Hess1, obj1->PObjData);
		obj2->DDFuncSqHess(Hess2, obj2->PObjData);
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
		HHess[0] = Hess1[0] + R[0][0] * RTH2[0] + R[1][0] * RTH2[1] + R[2][0] * RTH2[2];
		HHess[1] = Hess1[1] + R[0][1] * RTH2[0] + R[1][1] * RTH2[1] + R[2][1] * RTH2[2];
		HHess[2] = Hess1[2] + R[0][2] * RTH2[0] + R[1][2] * RTH2[1] + R[2][2] * RTH2[2];
		HHess[4] = Hess1[3] + R[0][1] * RTH2[3] + R[1][1] * RTH2[4] + R[2][1] * RTH2[5];
		HHess[5] = Hess1[4] + R[0][2] * RTH2[3] + R[1][2] * RTH2[4] + R[2][2] * RTH2[5];
		HHess[7] = Hess1[5] + R[0][2] * RTH2[6] + R[1][2] * RTH2[7] + R[2][2] * RTH2[8];
		double H2v[3];
		H2v[0] = Hess2[0] * v[0] + Hess2[1] * v[1] + Hess2[2] * v[2];
		H2v[1] = Hess2[1] * v[0] + Hess2[3] * v[1] + Hess2[4] * v[2];
		H2v[2] = Hess2[2] * v[0] + Hess2[4] * v[1] + Hess2[5] * v[2];

		HHess[3] = 0, HHess[6] = 0, HHess[8] = 0;
		for (int j = 0;j < 3;j++) {
			HHess[3] += R[j][0] * H2v[j];
			HHess[6] += R[j][1] * H2v[j];
			HHess[8] += R[j][2] * H2v[j];
		}
		HHess[9] = v[0] * H2v[0] + v[1] * H2v[1] + v[2] * H2v[2];
		/*-------------------------HessOfVCCD(obj1, obj2, x, HHess, R, p0, v, t);------------------------*/


		if (!CholeskySolveNewtonStepMove(HHess, HGrad, q)) {		//Solve the Newton step p=-Hess^{+}Grad
			if (!EigenDecompSolveNewtonStepMove(HHess, HGrad, q)) {
				Result->if_success = 0;
				return 0;
			}
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

		//double V_x;
		//ValueOfVCCD(obj1, obj2, x, &V_x, R, p0, v, t);
		/*-------------------------ValueOfVCCD(obj1, obj2, x, &V_x, R, p0, v, t);------------------------*/
		/*
		double Val1;
		obj1->DDFuncSqValu(&Val1, obj1->PObjData);
		double Val2;
		obj2->DDFuncSqValu(&Val2, obj2->PObjData);
		V_x = Val1 + Val2;*/
		/*-------------------------ValueOfVCCD(obj1, obj2, x, &V_x, R, p0, v, t);------------------------*/
		
		double q_in_2[3] = { 0,0,0 };
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				q_in_2[i] += R[i][j] * q[j];
			}
		}

		double alpha = 1, ro = 0.8;
		double V_al;
		double x_al[3], x_in_2_al[3], t_al;
		t_al = t + q[3];
		x_al[0] = x[0] + q[0];
		x_al[1] = x[1] + q[1];
		x_al[2] = x[2] + q[2];		
		x_in_2_al[0] = x_in_2[0] + q_in_2[0] - (t_al - t) * v2_in_2[0];
		x_in_2_al[1] = x_in_2[1] + q_in_2[1] - (t_al - t) * v2_in_2[1];
		x_in_2_al[2] = x_in_2[2] + q_in_2[2] - (t_al - t) * v2_in_2[2];
		//ValueOfVatXMove(obj1, obj2, x_al, &V_al, R, p0, v, t_al);
		obj1->DDFuncSqValuIndep(x_al, &Val1, obj1->PObjData);
		obj2->DDFuncSqValuIndep(x_in_2_al, &Val2, obj2->PObjData);
		V_al = Val1 + Val2;
		
		double  alpha_b = ro;
		double V_al_b;
		double x_al_b[3], x_in_2_al_b[3], t_al_b;
		t_al_b = t + alpha_b * q[3];
		x_al_b[0] = x[0] + alpha_b * q[0];
		x_al_b[1] = x[1] + alpha_b * q[1];
		x_al_b[2] = x[2] + alpha_b * q[2];
		x_in_2_al_b[0] = x_in_2[0] + alpha_b * q_in_2[0] - (t_al_b - t) * v2_in_2[0];
		x_in_2_al_b[1] = x_in_2[1] + alpha_b * q_in_2[1] - (t_al_b - t) * v2_in_2[1];
		x_in_2_al_b[2] = x_in_2[2] + alpha_b * q_in_2[2] - (t_al_b - t) * v2_in_2[2];
		//ValueOfVatXMove(obj1, obj2, x_al_b, &V_al_b, R, p0, v, t_al_b);
		obj1->DDFuncSqValuIndep(x_al_b, &Val1, obj1->PObjData);
		obj2->DDFuncSqValuIndep(x_in_2_al_b, &Val2, obj2->PObjData);
		V_al_b = Val1 + Val2;

		while (V_al - V_al_b > 0 && alpha > 0.01) {
			alpha = alpha_b;
			V_al = V_al_b;
			alpha_b *= ro;
			
			t_al = t_al_b;
			x_al[0] = x_al_b[0];
			x_al[1] = x_al_b[1];
			x_al[2] = x_al_b[2];
			x_in_2_al[0] = x_in_2_al_b[0];
			x_in_2_al[1] = x_in_2_al_b[1];
			x_in_2_al[2] = x_in_2_al_b[2];

			t_al_b = t + alpha_b * q[3];
			x_al_b[0] = x[0] + alpha_b * q[0];
			x_al_b[1] = x[1] + alpha_b * q[1];
			x_al_b[2] = x[2] + alpha_b * q[2];
			x_in_2_al_b[0] = x_in_2[0] + alpha_b * q_in_2[0] - (t_al_b - t) * v2_in_2[0];
			x_in_2_al_b[1] = x_in_2[1] + alpha_b * q_in_2[1] - (t_al_b - t) * v2_in_2[1];
			x_in_2_al_b[2] = x_in_2[2] + alpha_b * q_in_2[2] - (t_al_b - t) * v2_in_2[2];
			//ValueOfVatXMove(obj1, obj2, x_al_b, &V_al_b, R, p0, v, t_al_b);
			obj1->DDFuncSqValuIndep(x_al_b, &Val1, obj1->PObjData);
			obj2->DDFuncSqValuIndep(x_in_2_al_b, &Val2, obj2->PObjData);
			V_al_b = Val1 + Val2;
		}
		t = t_al;
		x[0] = x_al[0];
		x[1] = x_al[1];
		x[2] = x_al[2]; 
		x_in_2[0] = x_in_2_al[0];
		x_in_2[1] = x_in_2_al[1];
		x_in_2[2] = x_in_2_al[2];
		/*
		for (int i = 0;i < 3;i++) {
			p[i] = p0[i] + t * v[i];
		}
		x_in_2[0] = 0, x_in_2[1] = 0, x_in_2[2] = 0;
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				x_in_2[i] += R[i][j] * x[j];
			}
			x_in_2[i] += p[i];
		}
		*/
		//GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);
		/*-------------------------GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);------------------------*/
		obj1->DDFuncSqGrad(x, Grad1, obj1->PObjData);
		obj2->DDFuncSqGrad(x_in_2, Grad2, obj2->PObjData);

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
		/*-------------------------GradOfVCCD(obj1, obj2, x, HGrad, R, p0, v, t);------------------------*/
		GradNorm = sqrt(HGrad[0] * HGrad[0] + HGrad[1] * HGrad[1] + HGrad[2] * HGrad[2] + HGrad[3] * HGrad[3]);

		icount++;
	}
	if (icount > 49) {
		Result->if_success = 0;
		/*
		cout << "--------------------------------------------------" << endl;
		cout << "Failed: icount = " << icount << " > 49" << endl;
		cout << "Failed V: --------------" << endl;
		double V_fail;
		obj1->DDFuncSqValu(&Val1, obj1->PObjData);
		obj2->DDFuncSqValu(&Val2, obj2->PObjData);
		V_fail = Val1 + Val2;
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
		*/
	}
	else {
		Result->if_success = 1;
	}
	obj1->DDFuncSqValu(&Val1, obj1->PObjData);
	obj2->DDFuncSqValu(&Val2, obj2->PObjData);
	Result->V_min = Val1 + Val2;
	//ValueOfVatXMove(obj1, obj2, x, &(Result->V_min), R, p0, v, t);
	if (Result->V_min < ERROR_COLLIDE_V) {
		Result->if_collide = 1;
	}
	else {
		Result->if_collide = 0;
	}
	Result->num_iter = icount;

	return 1;
}