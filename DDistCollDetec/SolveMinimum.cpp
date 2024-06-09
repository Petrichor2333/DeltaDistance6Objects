#include "pch.h"
#include "config.h"
#include "SolvingMinimum.h"
#include "DDFunc.h"

#undef abs
/*
void ValueOfVatX(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, double* x, double* val, double R[][3], double* p) {
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
}*/


bool CholeskySolveNewtonStep(double* H, double* g, double* x) {
	double L[6];
	double l0_sq, l3_sq, l5_sq;

	l0_sq = H[0];
	if (l0_sq < fabs(l0_sq)*epsilon) { return 0; }
	L[0] = sqrt(l0_sq);
	L[1] = H[1] / L[0];
	L[2] = H[2] / L[0];
	l3_sq = H[3] - L[1] * L[1];
	if (l3_sq < fabs(l3_sq)*epsilon) { return 0; }
	L[3] = sqrt(l3_sq);
	L[4] = (H[4] - L[1] * L[2]) / L[3];
	l5_sq = H[5] - L[2] * L[2] - L[4] * L[4];
	if (l5_sq < fabs(l5_sq)*epsilon) { return 0; }
	L[5] = sqrt(l5_sq);

	double y[3];
	y[0] = -g[0] / L[0];										//x=-Hess^{-1}Grad, so "-g"
	y[1] = (-g[1] - L[1] * y[0]) / L[3];
	y[2] = (-g[2] - L[2] * y[0] - L[4] * y[1]) / L[5];

	x[2] = y[2] / L[5];								
	x[1] = (y[1] - L[4] * x[2]) / L[3];
	x[0] = (y[0] - L[1] * x[1] - L[2] * x[2]) / L[0];
	return 1;
}

bool EigenDecompSolveNewtonStep(double* H, double* g, double* x) {
	integer mat_size = 3;
	char JOBZ = 'V';
	char UPLO = 'U';
	double P[9] = { H[0], H[1], H[2],
					H[1], H[3], H[4],
					H[2], H[4], H[5] };
	integer LDA = 3;
	double eig_val[3];
	double WORK[9];
	integer LWORK = 9;
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
	double y[3], z[3];
	for (int i = 0;i < 3;i++) {
		y[i] = 0;
		for (int j = 0;j < 3;j++) {
			y[i] -= P[3 * i + j] * g[j];
		}
	}

	//	double tolerance = epsilon * fabs(eig_val[2]);
	for (int i = 0;i < 3;i++) {
		z[i] = eig_val[i] > epsilon * fabs(eig_val[i]) ? y[i] / fabs(eig_val[i]) : 0;
	}

	for (int i = 0;i < 3;i++) {
		x[i] = 0;
		for (int j = 0;j < 3;j++) {
			x[i] += P[3 * j + i] * z[j];
		}
	}
	return 1;
}

bool SolveForCollision(ObjectDDFuncInfo* obj1, ObjectDDFuncInfo* obj2, ResultsOfSolvingMinimum * Result, double R1[][3], double* p1, double R2[][3], double* p2) {
	double R[3][3];
	double p[3];
	double x[3];

	for (int i = 0;i < 3;i++) {
		p[i] = 0;
		for (int j = 0;j < 3;j++) {
			p[i] += R2[j][i] * (p1[j] - p2[j]);
			R[i][j] = 0;
			for (int k = 0;k < 3;k++) {
				R[i][j] += R2[k][i] * R1[k][j];
			}
		}
	}

	for (int i = 0;i < 3;i++) {
		x[i] = 0;
		for (int j = 0;j < 3;j++) {
			x[i] -= R[j][i] * p[j];
		}
		x[i] *= 0.5;
	}

	double Val1, Val2;

	double Grad[3];
	double Grad1[3], Grad2[3];	
	double x_in_2[3];
	for (int i = 0;i < 3;i++) {
		x_in_2[i] = p[i];
		for (int j = 0;j < 3;j++) {
			x_in_2[i] += R[i][j] * x[j];
		}		
	}	
	//GradOfV(obj1, obj2, x, Grad, R, p);
	/*----------------GradOfV(obj1, obj2, x, Grad, R, p);----------------*/
	obj1->DDFuncSqGrad(x, Grad1, obj1->PObjData);
	obj2->DDFuncSqGrad(x_in_2, Grad2, obj2->PObjData);
	for (int i = 0;i < 3;i++) {
		Grad[i] = Grad1[i];
		for (int j = 0;j < 3;j++) {
			Grad[i] += R[j][i] * Grad2[j];
		}
	}
	/*----------------GradOfV(obj1, obj2, x, Grad, R, p);----------------*/

	double GradNorm = sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]);
	
	double q[3];	//Newton Step

	int icount = 0;

	while (GradNorm > ERROR_GRAD && icount < 51)  {	
		double Hess[6];
		//HessOfV(obj1, obj2, x, Hess, R, p);
		/*----------------HessOfV(obj1, obj2, x, Hess, R, p);----------------*/
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
		Hess[0] = Hess1[0] + R[0][0] * RTH2[0] + R[1][0] * RTH2[1] + R[2][0] * RTH2[2];
		Hess[1] = Hess1[1] + R[0][1] * RTH2[0] + R[1][1] * RTH2[1] + R[2][1] * RTH2[2];
		Hess[2] = Hess1[2] + R[0][2] * RTH2[0] + R[1][2] * RTH2[1] + R[2][2] * RTH2[2];
		Hess[3] = Hess1[3] + R[0][1] * RTH2[3] + R[1][1] * RTH2[4] + R[2][1] * RTH2[5];
		Hess[4] = Hess1[4] + R[0][2] * RTH2[3] + R[1][2] * RTH2[4] + R[2][2] * RTH2[5];
		Hess[5] = Hess1[5] + R[0][2] * RTH2[6] + R[1][2] * RTH2[7] + R[2][2] * RTH2[8];
		/*----------------HessOfV(obj1, obj2, x, Hess, R, p);----------------*/

		if (!CholeskySolveNewtonStep(Hess, Grad, q)) {		//Solve the Newton step p = -Hess^{+}Grad
			if (!EigenDecompSolveNewtonStep(Hess, Grad, q)) {
				Result->if_success = 0;
				return 0;
			}
		}
		/*
		if (!SolveNewtonStep(Hess, Grad, q)) {
			Result->if_success = 0; 
			return 0; 
		}*/	
		
		double ModPara = 0.1;
		
		double q_norm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		double max_step = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])/2;
		if (q_norm > max_step) {
			q[0] *= (max_step / q_norm);
			q[1] *= (max_step / q_norm);
			q[2] *= (max_step / q_norm);
			q_norm = max_step;
		}		

		//double V_x;
		//ValueOfV(obj1, obj2, x, &V_x, R, p);
		/*----------------ValueOfV(obj1, obj2, x, &V_x, R, p);----------------*/
		/*
		double Val1;
		obj1->DDFuncSqValu(&Val1, obj1->PObjData);
		double Val2;
		obj2->DDFuncSqValu(&Val2, obj2->PObjData);
		V_x = Val1 + Val2;*/
		/*----------------ValueOfV(obj1, obj2, x, &V_x, R, p);----------------*/
		
		double q_in_2[3] = { 0,0,0 };
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				q_in_2[i] += R[i][j] * q[j];
			}
		}
		
		double alpha = 1, ro = 0.7;
		double V_al;
		double x_al[3], x_in_2_al[3];
		x_al[0] = x[0] + q[0];
		x_al[1] = x[1] + q[1];
		x_al[2] = x[2] + q[2];
		x_in_2_al[0] = x_in_2[0] + q_in_2[0];
		x_in_2_al[1] = x_in_2[1] + q_in_2[1];
		x_in_2_al[2] = x_in_2[2] + q_in_2[2];
		//ValueOfVatX(obj1, obj2, x_al, &V_al, R, p);
		obj1->DDFuncSqValuIndep(x_al, &Val1, obj1->PObjData);
		obj2->DDFuncSqValuIndep(x_in_2_al, &Val2, obj2->PObjData);
		V_al = Val1 + Val2;

		double alpha_b = ro;
		double V_al_b;		
		double x_al_b[3], x_in_2_al_b[3];
		x_al_b[0] = x[0] + alpha_b * q[0];
		x_al_b[1] = x[1] + alpha_b * q[1];
		x_al_b[2] = x[2] + alpha_b * q[2];
		x_in_2_al_b[0] = x_in_2[0] + alpha_b * q_in_2[0];
		x_in_2_al_b[1] = x_in_2[1] + alpha_b * q_in_2[1];
		x_in_2_al_b[2] = x_in_2[2] + alpha_b * q_in_2[2];
		//ValueOfVatX(obj1, obj2, x_al_b, &V_al_b, R, p);
		obj1->DDFuncSqValuIndep(x_al_b, &Val1, obj1->PObjData);
		obj2->DDFuncSqValuIndep(x_in_2_al_b, &Val2, obj2->PObjData);
		V_al_b = Val1 + Val2;

		while ( V_al > V_al_b && alpha > 0.01) {
			alpha = alpha_b;
			V_al = V_al_b;
			alpha_b *= ro;

			x_al[0] = x_al_b[0];
			x_al[1] = x_al_b[1];
			x_al[2] = x_al_b[2];
			x_in_2_al[0] = x_in_2_al_b[0];
			x_in_2_al[1] = x_in_2_al_b[1];
			x_in_2_al[2] = x_in_2_al_b[2];	
			
			x_al_b[0] = x[0] + alpha_b * q[0];
			x_al_b[1] = x[1] + alpha_b * q[1];
			x_al_b[2] = x[2] + alpha_b * q[2];
			x_in_2_al_b[0] = x_in_2[0] + alpha_b * q_in_2[0];
			x_in_2_al_b[1] = x_in_2[1] + alpha_b * q_in_2[1];
			x_in_2_al_b[2] = x_in_2[2] + alpha_b * q_in_2[2];
			//ValueOfVatX(obj1, obj2, x_al_b, &V_al_b, R, p);
			obj1->DDFuncSqValuIndep(x_al_b, &Val1, obj1->PObjData);
			obj2->DDFuncSqValuIndep(x_in_2_al_b, &Val2, obj2->PObjData);
			V_al_b = Val1 + Val2;
			
		}
		x[0] = x_al[0];
		x[1] = x_al[1];
		x[2] = x_al[2];
		x_in_2[0] = x_in_2_al[0];
		x_in_2[1] = x_in_2_al[1];
		x_in_2[2] = x_in_2_al[2];
		/*
		x_in_2[0] = 0, x_in_2[1] = 0, x_in_2[2] = 0;
		for (int i = 0;i < 3;i++) {
			for (int j = 0;j < 3;j++) {
				x_in_2[i] += R[i][j] * x[j];
			}
			x_in_2[i] += p[i];
		}*/

		//GradOfV(obj1, obj2, x, Grad, R, p);
		/*----------------GradOfV(obj1, obj2, x, Grad, R, p);----------------*/
		obj1->DDFuncSqGrad(x, Grad1, obj1->PObjData);
		obj2->DDFuncSqGrad(x_in_2, Grad2, obj2->PObjData);
		for (int i = 0;i < 3;i++) {
			Grad[i] = Grad1[i];
			for (int j = 0;j < 3;j++) {
				Grad[i] += R[j][i] * Grad2[j];
			}
		}
		/*----------------GradOfV(obj1, obj2, x, Grad, R, p);----------------*/

		GradNorm = sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]);

		icount++;
	}
	
	if (icount > 49) {
		Result->if_success = 0;
		/*
		cout << "--------------------------------------------------" << endl;
		cout << "Failed: icount = "<< icount<<" > 49" << endl;
		cout << "Failed V: --------------" << endl;
		double V_fail;
		obj1->DDFuncSqValu(&Val1, obj1->PObjData);
		obj2->DDFuncSqValu(&Val2, obj2->PObjData);
		V_fail = Val1 + Val2;
		//ValueOfVatX(obj1, obj2, x, &V_fail, R, p);
		cout << V_fail << endl;
		cout << endl;
		cout << "Failed Grad:" << endl;
		cout << Grad[0] << ", " << Grad[1] << ", " << Grad[2] << endl;
		cout << endl;
		cout << "Failed p:" << endl;
		cout << q[0] << ", " << q[1] << ", " << q[2] << endl;
		cout << endl;
		cout << "x:" << endl;
		cout << x[0] << ", " << x[1] << ", " << x[2] << endl;
		cout << "--------------------------------------------------" << endl;
		*/
	}
	else {
		Result->if_success = 1;
	}	
	obj1->DDFuncSqValu(&Val1, obj1->PObjData);
	obj2->DDFuncSqValu(&Val2, obj2->PObjData);
	Result->V_min = Val1 + Val2;
	//ValueOfVatX(obj1, obj2, x, &(Result->V_min), R, p);
	Result->if_collide = Result->V_min > ERROR_COLLIDE_V ? 1 : 0;
	Result->num_iter = icount;
	Result->ConvPoint[0] = x[0], Result->ConvPoint[1] = x[1], Result->ConvPoint[2] = x[2];
	return 1;
}

